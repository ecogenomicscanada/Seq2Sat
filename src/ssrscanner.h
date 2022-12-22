#ifndef SSRSCANNER_H
#define SSRSCANNER_H

#include <string>
#include <algorithm>
#include <vector>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <bits/stdc++.h>

#include "options.h"
#include "read.h"
#include "genotype.h"
#include "util.h"
#include "edlib.h"
#include "common.h"
#include "editdistance.h"

using namespace std;

class SsrScanner {
public:
    SsrScanner(Options* opt);
    ~SsrScanner();
    
    std::string scanVar(Read* & r1);
    bool scanVar(Read* & r1, Read* & r2);
    
    static std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> report(Options * & mOptions, std::map<std::string, std::map<std::string, Genotype>> & allGenotypeMap);
    static std::map<std::string, std::map<std::string, Genotype>> merge(std::vector<std::map<std::string, std::map<std::string, Genotype>>> & totalGenotypeSsrMapVec);
    void static merge(std::vector<std::map<std::string, std::map<std::string, int>>> & totalSexLocVec, Options * & mOptions);
    inline std::map<std::string, std::map< std::string, Genotype>> getGenotypeMap() {return tmpAllGenotypeMap;};
    static std::map<int, std::string> doSimpleAlignment(Options * & mOptions, const char* & qData, int qLength, const char* & tData, int tLength);
    //inline std::map<std::string, std::vector<std::pair<std::string, Genotype>>> getSortedGenotypeMap() {return sortedAllGenotypeMap;};
    inline std::map<std::string, std::map<std::string, int>> getSexLoc(){return tmpSexMap;};
    
private:
    
    int doAlignment(const char* & qData, int qLength, const std::string & qName,
                     const char* & tData, int tLength, const std::string & tName, 
                     Variance & variance, bool printAlignment);

    void doScanVariance(EdlibAlignResult & result, Variance & variance, const char* & qData, const char* & tData, const int position);
    
    void printVariance(EdlibAlignResult & result, Variance & variance, 
                       const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position);
    
    int doSubStrAlignment(const char* & pData, int pLength, const char* & rData, int rLength,
                          int & numMismaches, int & numDeletions, int & numInsertions, bool rev);
    
    void printAlignment(const char* & query, const char* & target,
            const unsigned char*  alignment, const int alignmentLength,
            const int position, const EdlibAlignMode modeCode);
    
    std::pair<size_t, int> analyzeMRA(std::string rawStr,  const std::string & ssr, std::size_t & ffTrimPos, std::size_t & rfTrimPos);//rawStr must be passed by value not reference
    
    void preAnalyze(Read* & r1, std::size_t & ffpos, std::size_t & rfpos, bool & mraAnalyze);
    
    void analyzeMRASub(Genotype & genotype, std::map<int, std::string> & subMap);
    
    void analyzeFlankingRegion(Genotype & genotype, 
                               Variance & ffVar, int endPosF, std::pair<size_t, int> & mraPosLenReadF,
                               Variance & rfVar, int endPosR, std::pair<size_t, int> & mraPosLenReadR);
    
    void analyzeUnmatchedSeq(std::string & unmatchedStr, bool ff);
    
    bool analyzeMRAVar(Genotype & genotype, Variance & fVar, Variance & rVar);
    
    bool checkMRASub(Variance & fVar);
    
    bool checkMRAIndel(Variance & fVar);
    
    void resetData();
    
    //std::size_t mutationMatch(std::string target, std::string query, bool rev = false);//must passed by value for target;
    
    std::size_t mutationMatchFR(std::string target, std::string query, int mis, bool rev = false);
    
    //std::string getGenotype(std::string & mra, std::string & ssr);
    
    std::string trim2ends(Read* & tmpRead, int startPos, int endPos);
    
    std::vector<std::pair<std::string, Genotype>> sortGenotypeMap(std::map<std::string, Genotype> & genoMap);
    

private:
    Options* mOptions;
    const char* fpData;
    int fpLength;
    const char* rpData;
    int rpLength;
    const char* target;
    int targetLength;
    const char* readSeq;
    int readLength;
    priority_queue<int> bestScores;
    LocVar* locVarIt;
    std::string readName;
    std::stringstream ss;
    std::map<std::string, std::map<std::string, Genotype>> tmpAllGenotypeMap;//marker, seq, geno
    //Sex tmpSex;
    std::map<std::string, std::map<std::string, int>> tmpSexMap;
    //std::map<std::string, std::vector<std::pair<std::string, Genotype>>> sortedAllGenotypeMap;
    int ffEndPos;
    int rfStartPos;
    std::string readWithoutFF;
    std::string readWithoutFFRF;
    std::string locMra;
    std::string enhancer;
    std::string returnedlocus;
    bool checkLoci;
    uint32 mismachesFF;
    uint32 mismachesRF;
    int minReadLen;
};

#endif /* SSRSCANNER_H */

