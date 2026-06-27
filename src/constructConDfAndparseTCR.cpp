// Rcpp replacement for .parseTCR
// By Qile Yang

#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "scRepHelper.h"

#define BarcodeIndciesMap std::unordered_map<std::string, std::vector<int>>

std::vector<std::vector<int>> constructBarcodeIndex(
    std::vector<std::string>& conDfBarcodes, std::vector<std::string> data2Barcodes
) {
    std::vector<std::vector<int>> outputBarcodeIndex (conDfBarcodes.size());
    BarcodeIndciesMap data2BarcodeIndiciesMap = scRepHelper::stringToIndiciesMap(data2Barcodes);

    for (int i = 0; i < (int) conDfBarcodes.size(); i++) {
        std::string& barcode = conDfBarcodes[i];
        if (data2BarcodeIndiciesMap.find(barcode) != data2BarcodeIndiciesMap.end()) {
            outputBarcodeIndex[i] = data2BarcodeIndiciesMap[barcode];
        }
    }
    return outputBarcodeIndex;
}

class TcrParser {
public:
    // variable for the eventual output Con.df
    std::vector<std::vector<std::string>> conDf;

    // variables for *references* to columns on data2
    Rcpp::CharacterVector data2ChainTypes;
    Rcpp::CharacterVector data2Tcr1;
    Rcpp::CharacterVector data2Tcr2;
    Rcpp::CharacterVector data2Cdr3;
    Rcpp::CharacterVector data2Cdr3Nt;

    // variable for helper barcode index
    std::vector<std::vector<int>> barcodeIndex;

    // optional full-length sequence columns carried in the same pass.
    // For each retained data2 column we keep one output slot per chain
    // (slot 1 = TRA/TRG, slot 2 = TRB/TRD), ';'-joined to mirror the
    // multi-contig concatenation of cdr3_nt1 / cdr3_nt2.
    std::vector<std::string> seqColNames;
    std::vector<Rcpp::CharacterVector> seqCols;
    std::vector<std::vector<std::string>> seqSlot1;
    std::vector<std::vector<std::string>> seqSlot2;

    // constructor
    TcrParser(
        Rcpp::DataFrame& data2, std::vector<std::string>& uniqueData2Barcodes,
        std::vector<std::string> retainCols
    ) {
        // construct conDf, initializing the matrix to "NA" *strings*
        conDf = scRepHelper::initStringMatrix(
            7, uniqueData2Barcodes.size(), "NA"
        );
        conDf[0] = uniqueData2Barcodes;

        // set references to fixed data2 columns
        data2ChainTypes = data2[data2.findName("chain")];
        data2Cdr3 = data2[data2.findName("cdr3")];
        data2Cdr3Nt = data2[data2.findName("cdr3_nt")];

        // setting reference to the TCR columns assuming all extra columns come before
        data2Tcr1 = data2[data2.findName("TCR1")];
        data2Tcr2 = data2[data2.findName("TCR2")];

        // construct barcodeIndex
        barcodeIndex = constructBarcodeIndex(
            uniqueData2Barcodes, data2[data2.findName("barcode")]
        );

        // set up optional retained sequence columns
        seqColNames = retainCols;
        for (int c = 0; c < (int) seqColNames.size(); c++) {
            seqCols.push_back(data2[data2.findName(seqColNames[c])]);
            seqSlot1.push_back(std::vector<std::string>(uniqueData2Barcodes.size(), "NA"));
            seqSlot2.push_back(std::vector<std::string>(uniqueData2Barcodes.size(), "NA"));
        }
    }

    // Rcpp implementation of .parseTCR()
    TcrParser& parseTCR() {
        for (int y = 0; y < (int) conDf[0].size(); y++) {
            for (int index : barcodeIndex[y]) {
                std::string chainType = std::string(data2ChainTypes[index]);
                if (chainType == "TRA" || chainType == "TRG") {
                    handleTcr1(y, index);
                } else if (chainType == "TRB" || chainType == "TRD") {
                    handleTcr2(y, index);
                } else {
                    Rcpp::stop("Invalid chain type: " + chainType + " for barcode: " + conDf[0][y]);
                }
            }
        }
        return *this;
    }

    // parseTCR() helpers

    void handleTcr1(int y, int data2index) {
        handleTcr(y, data2index, data2Tcr1, 1, 2, 3, 1);
    }

    void handleTcr2(int y, int data2index) {
        handleTcr(y, data2index, data2Tcr2, 4, 5, 6, 2);
    }

    void handleTcr(
        int y, int data2index, Rcpp::CharacterVector& data2tcr,
        int tcr, int cdr3aa, int cdr3nt, int slot
    ) {
        // Whether this is the first contig for this chain slot determines
        // set-vs-append; the retained sequence slots must use the SAME
        // decision so their ';'-join order matches cdr3_nt exactly.
        bool firstChain = (conDf[tcr][y] == "NA");
        if (firstChain) {
            conDf[tcr][y] = data2tcr[data2index];
            conDf[cdr3aa][y] = data2Cdr3[data2index];
            conDf[cdr3nt][y] = data2Cdr3Nt[data2index];
        } else {
            conDf[tcr][y] += ";" + data2tcr[data2index];
            conDf[cdr3aa][y] += ";" + data2Cdr3[data2index];
            conDf[cdr3nt][y] += ";" + data2Cdr3Nt[data2index];
        }

        for (int c = 0; c < (int) seqCols.size(); c++) {
            std::vector<std::string>& slotVec = (slot == 1) ? seqSlot1[c] : seqSlot2[c];
            std::string val = std::string(seqCols[c][data2index]);
            if (firstChain) {
                slotVec[y] = val;
            } else {
                slotVec[y] += ";" + val;
            }
        }
    }

    // return Con.df after TCR parsing
    Rcpp::DataFrame getConDf() {
        if (seqColNames.empty()) {
            return Rcpp::DataFrame::create(
                Rcpp::Named("barcode") = conDf[0],
                Rcpp::Named("TCR1") = conDf[1],
                Rcpp::Named("cdr3_aa1") = conDf[2],
                Rcpp::Named("cdr3_nt1") = conDf[3],
                Rcpp::Named("TCR2") = conDf[4],
                Rcpp::Named("cdr3_aa2") = conDf[5],
                Rcpp::Named("cdr3_nt2") = conDf[6]
            );
        }

        // Dynamic column set: 7 base columns + 2 per retained sequence column.
        int nBase = 7;
        int nExtra = 2 * (int) seqColNames.size();
        Rcpp::List cols(nBase + nExtra);
        Rcpp::CharacterVector names(nBase + nExtra);

        const char* baseNames[7] = {
            "barcode", "TCR1", "cdr3_aa1", "cdr3_nt1", "TCR2", "cdr3_aa2", "cdr3_nt2"
        };
        for (int i = 0; i < nBase; i++) {
            cols[i] = Rcpp::wrap(conDf[i]);
            names[i] = baseNames[i];
        }
        for (int c = 0; c < (int) seqColNames.size(); c++) {
            int i1 = nBase + 2 * c;
            int i2 = i1 + 1;
            cols[i1] = Rcpp::wrap(seqSlot1[c]);
            names[i1] = seqColNames[c] + "1";
            cols[i2] = Rcpp::wrap(seqSlot2[c]);
            names[i2] = seqColNames[c] + "2";
        }

        cols.attr("names") = names;
        cols.attr("row.names") = Rcpp::IntegerVector::create(
            NA_INTEGER, -(int) conDf[0].size()
        );
        cols.attr("class") = "data.frame";
        return Rcpp::as<Rcpp::DataFrame>(cols);
    }
};

// [[Rcpp::export]]
Rcpp::DataFrame rcppConstructConDfAndParseTCR(
    Rcpp::DataFrame& data2, std::vector<std::string> uniqueData2Barcodes
) {
    std::vector<std::string> noCols;
    return TcrParser(data2, uniqueData2Barcodes, noCols).parseTCR().getConDf();
}

// [[Rcpp::export]]
Rcpp::DataFrame rcppConstructConDfAndParseTCRWithSeqs(
    Rcpp::DataFrame& data2, std::vector<std::string> uniqueData2Barcodes,
    std::vector<std::string> retainCols
) {
    return TcrParser(data2, uniqueData2Barcodes, retainCols).parseTCR().getConDf();
}
