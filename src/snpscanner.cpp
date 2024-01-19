#include <valarray>

#include "snpscanner.h"

SnpScanner::SnpScanner(Options* opt) {
    mOptions = opt;
    //subGenotypeMap.clear();
    subSeqsMap.clear();
    target = NULL;
    targetLength = 0;
    readSeq = NULL;
    readLength = 0;
    readName = "";
    fpData = NULL;
    fpLength = 0;
    rpData = NULL;
    rpLength = 0;
    ss.str();
}

SnpScanner::SnpScanner(const SnpScanner& orig) {
}

SnpScanner::~SnpScanner() {
}

bool SnpScanner::scanVar(Read* & r1, Read* & r2) {
    return true;
}

std::string SnpScanner::scanVar(Read* & r1) {
    ss.str("");
    returnedlocus.clear();
    readSeq = r1->mSeq.mStr.c_str();
    readLength = r1->mSeq.mStr.length();
    readName = r1->mName;
    std::map<std::string, MatchTrim> locMap;
    
    for (auto & it : mOptions->mLocSnps.refLocMap) {
        
        if(r1->mSeq.length() < (it.second.fp.length() + it.second.ft.length() + it.second.ref.length() + it.second.rt.length() + it.second.rp.length())){
            continue;
        }
        
        bool goRP = false;
        int trimF = 0;
        int fpMismatches = (int)edit_distance(it.second.fp.mStr, r1->mSeq.mStr.substr(0, it.second.fp.length()));
        if (fpMismatches <= mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq) {
            trimF = it.second.fp.length();
            goRP = true;
        } else {
            fpData = it.second.fp.mStr.c_str();
            fpLength = it.second.fp.length();
            auto endBoolF = doPrimerAlignment(fpData, fpLength, mOptions->mSex.sexMarker, readSeq, readLength, r1->mName, true);
             if (get<2>(endBoolF) && (get<1>(endBoolF) <= r1->length())) {
                 fpMismatches = get<0>(endBoolF);//already checked for mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq
                 if((get<1>(endBoolF) + it.second.ft.length() + it.second.ref.length() + it.second.rt.length() + it.second.rp.length()) <= r1->mSeq.length()){
                     trimF = get<1>(endBoolF);
                     goRP = true;
                 } else {
                     goRP = false;
                 }
             } else {
                 goRP = false;
             }
        }

        if (goRP) {
            MatchTrim mTrim;
            int rpMismatches = (int) edit_distance(it.second.rp.mStr, r1->mSeq.mStr.substr(r1->mSeq.length() - it.second.rp.length()));
            if (rpMismatches <= mOptions->mLocSnps.mLocSnpOptions.maxMismatchesPSeq) {
                mTrim.totMismatches = fpMismatches + rpMismatches;
                mTrim.trimF = trimF;
                mTrim.trimedRefLenth = r1->mSeq.length() - trimF - it.second.rp.mStr.length();
                locMap[it.second.name] = mTrim;
            } else {
                rpData = it.second.rp.mStr.c_str();
                rpLength = it.second.rp.length();
                auto endBoolR = doPrimerAlignment(rpData, rpLength, it.second.name, readSeq, readLength, r1->mName, true);
                if (get<2>(endBoolR) && (get<1>(endBoolR) <= r1->mSeq.mStr.length()) && 
                        ((trimF + it.second.ft.length() + it.second.ref.length() + it.second.rt.length() + it.second.rp.length()) <= get<1>(endBoolR))) {
                    mTrim.totMismatches = fpMismatches + rpMismatches;
                    mTrim.trimF = trimF;
                    mTrim.trimedRefLenth = get<1>(endBoolR) - trimF - it.second.rp.mStr.length();
                    locMap[it.second.name] = mTrim;
                }
            }
        }
    }

    if (locMap.empty()) {
        return returnedlocus;
    }

    std::string locName = "";
    if (locMap.size() == 1) {
        locName = locMap.begin()->first;
        //if(mOptions->debug) cCout("single value: " + locName, 'r');
    } else {
        std::vector<int> seqScoreVec;
        for (const auto & it : locMap) {
            seqScoreVec.push_back(it.second.totMismatches);
        }
        auto minValue = *std::min_element(seqScoreVec.begin(), seqScoreVec.end());
        //warning, what if there are multiple identical values
        seqScoreVec.clear();

        for (const auto & it : locMap) {
            if (it.second.totMismatches == minValue) {
                locName = it.first;
                ss.str();
                break; //warning, what if there are multiple identical values
            }
        }
    }
    
    locSnpIt = &(mOptions->mLocSnps.refLocMap[locName]);

    if (r1->mSeq.length() < (locSnpIt->fp.length() + locSnpIt->ft.length() + locSnpIt->ref.length() + locSnpIt->rt.length() + locSnpIt->rp.length())) {
        return returnedlocus;
    } else {
        if(locMap[locName].trimF < r1->mSeq.length()){
            r1->trimFront(locMap[locName].trimF);
            if(locMap[locName].trimedRefLenth <= r1->mSeq.length()){
                r1->resize(locMap[locName].trimedRefLenth);
            } else {
                return returnedlocus;
            }
        } else {
            return returnedlocus;
        }
    }
    locMap.clear();

    if (mOptions->debug) cCout("detected marker: " + locSnpIt->name, 'r');
    subSeqsMap[locSnpIt->name][r1->mSeq.mStr]++;
    returnedlocus = locName;
    if (mOptions->mEdOptions.printRes) {
        cCout(ss.str(), 'g');
    }
    ss.str();
    return returnedlocus;
}

std::pair<bool, std::map<int, std::pair<Sequence, Sequence>>> SnpScanner::doAlignment(Options * & mOptions, std::string readName, const char* & qData, int qLength, std::string targetName, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_NW,
            mOptions->mEdOptions.alignTask,
            NULL, 0));

    std::pair<bool, std::map<int, std::pair<Sequence, Sequence>>> snpsMapPair;

    if (result.status == EDLIB_STATUS_OK) {
        if (mOptions->mEdOptions.printRes) {
            Variance variance;
            doScanVariance(mOptions, result, variance, qData, tData, *(result.endLocations));
            if (mOptions->mEdOptions.printRes) {
                //printVariance(result, variance, qData, readName, tData, targetName, *(result.endLocations));
            }

            variance.cleanVar();
        }
        bool snps = true;
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                std::string s1(1, tData[i]);
                std::string s2(1, qData[i]);
                snpsMapPair.second[i] = std::make_pair(Sequence(s1), Sequence(s2));
            } else if (cur == EDLIB_EDOP_INSERT) {
                snps = snp && false;
            } else if (cur == EDLIB_EDOP_DELETE) {
                snps = snp && false;
            }
        }

        snpsMapPair.first = snps;
    }
    edlibFreeAlignResult(result);
    return snpsMapPair;
}

std::pair<bool, std::set<int>> SnpScanner::doAlignment2(Options * & mOptions, std::string readName, const char* & qData, int qLength, std::string targetName, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_NW,
            mOptions->mEdOptions.alignTask,
            NULL, 0));

    std::pair<bool, std::set<int>> snpsSetPair;

    if (result.status == EDLIB_STATUS_OK) {
        if (mOptions->mEdOptions.printRes) {
            Variance variance;
            doScanVariance(mOptions, result, variance, qData, tData, *(result.endLocations));
            if (mOptions->mEdOptions.printRes) {
                //printVariance(result, variance, qData, readName, tData, targetName, *(result.endLocations));
            }

            variance.cleanVar();
        }
        bool snps = true;
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsSetPair.second.insert(i);
            } else if (cur == EDLIB_EDOP_INSERT) {
                snps = snp && false;
            } else if (cur == EDLIB_EDOP_DELETE) {
                snps = snp && false;
            }
        }
        snpsSetPair.first = snps;
    }
    edlibFreeAlignResult(result);
    return snpsSetPair;
}

void SnpScanner::doScanVariance(Options * & mOptions, EdlibAlignResult & result, Variance & variance,
        const char* & qData, const char* & tData, const int position) {

    int ti = -1, qi = -1;
    std::string inStr, deStr, snpStr;

    if (result.alignmentLength < 2) {
        //if (mOptions->mLocVars.locVarOptions.printRes) ss << "You must have at least 2 bp!\n";
        return;
    }

    if (mOptions->mEdOptions.modeCode == EDLIB_MODE_HW) {
        ti = position;
        for (int i = 0; i < result.alignmentLength; i++) {
            if (result.alignment[i] != EDLIB_EDOP_INSERT) {
                ti--;
            }
        }
    }

    auto pro = result.alignment[0];
    auto cur = pro;

    if (pro == EDLIB_EDOP_MATCH) {
        ti++;
        qi++;
    } else if (pro == EDLIB_EDOP_INSERT) {
        qi++;
        inStr.push_back(qData[qi]);
        if (result.alignmentLength == 1) {
            //get<1>(*vTuple).emplace_back(ti, inStr);
            variance.insMap[ti + 1] = inStr;
            inStr.clear();
        }
    } else if (pro == EDLIB_EDOP_DELETE) {
        ti++;
        deStr.push_back(tData[ti]);
        if (result.alignmentLength == 1) {
            //get<2>(*vTuple).emplace_back(ti, deStr);
            variance.delMap[ti] = deStr;
            deStr.clear();
        }
    } else if (pro == EDLIB_EDOP_MISMATCH) {
        ti++;
        qi++;
        snpStr.push_back(tData[ti]);
        snpStr.push_back('|');
        snpStr.push_back(qData[qi]);
        //get<0>(*vTuple).emplace_back(std::make_pair(ti, snpStr));
        variance.subMap[ti] = snpStr;
        snpStr.clear();
    }

    for (int i = 1; i < result.alignmentLength; i++) {
        pro = result.alignment[i - 1];
        cur = result.alignment[i];

        if (cur == EDLIB_EDOP_MATCH) {

            ti++;
            qi++;

            if (pro == EDLIB_EDOP_MATCH) {

            } else if (pro == EDLIB_EDOP_INSERT) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti - inStr.length()] = inStr;
                inStr.clear();
            } else if (pro == EDLIB_EDOP_DELETE) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            } else if (pro == EDLIB_EDOP_MISMATCH) {

            }
        } else if (cur == EDLIB_EDOP_INSERT) {

            qi++;
            inStr.push_back(qData[qi]);

            if (pro == EDLIB_EDOP_DELETE) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            }

            if (i == result.alignmentLength - 1) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti] = inStr;
                inStr.clear();
            }


        } else if (cur == EDLIB_EDOP_DELETE) {
            ti++;
            deStr.push_back(tData[ti]);
            if (pro == EDLIB_EDOP_INSERT) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti - inStr.length()] = inStr;
                inStr.clear();
            }

            if (i == result.alignmentLength - 1) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            }


        } else if (cur == EDLIB_EDOP_MISMATCH) {

            ti++;
            qi++;

            snpStr.push_back(tData[ti]);
            snpStr.push_back('|');
            snpStr.push_back(qData[qi]);
            //get<0>(*vTuple).emplace_back(std::make_pair(ti, snpStr));
            variance.subMap[ti] = snpStr;
            snpStr.clear();

            if (pro == EDLIB_EDOP_INSERT) {
                //get<1>(*vTuple).emplace_back(std::make_pair(ti, inStr));
                variance.insMap[ti - inStr.length()] = inStr;
                inStr.clear();
            } else if (pro == EDLIB_EDOP_DELETE) {
                //get<2>(*vTuple).emplace_back(std::make_pair(ti - deStr.length(), deStr));
                variance.delMap[ti - deStr.length()] = deStr;
                deStr.clear();
            } else {

            }
        }
    }
}

//void SnpScanner::printVariance(Options * & mOptions, EdlibAlignResult & result, Variance & variance,
//        const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position) {
//
//    if (mOptions->mLocVars.locVarOptions.printRes) {
//        ss << "End location: " << *result.endLocations << " -> alignmentLength: " << result.alignmentLength << " with alignment distance: " << result.editDistance << "\n";
//        ss << "Ref: " << tName << " -> Query: " << qName << "\n";
//    }
//    for (int i = 0; i < result.alignmentLength; i++) {
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << static_cast<unsigned> (result.alignment[i]);
//    }
//    if (mOptions->mLocVars.locVarOptions.printRes) ss << "\n";
//
//
//    for (const auto & it : variance.subMap) {
//        //std::string str = "snp: " + std::to_string(it.first) + " -> " + it.second;
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "snp: " << it.first << " -> " << it.second << "\n";
//    }
//
//    //get<0>(*vTuple).clear();
//    //get<0>(*vTuple).shrink_to_fit();
//
//    for (const auto & it : variance.insMap) {
//        //std::string str = "insertion: " + std::to_string(it.first) + " -> " + it.second;
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "insertion: " << it.first << " -> " << it.second << "\n";
//    }
//
//    //get<1>(*vTuple).clear();
//    //get<1>(*vTuple).shrink_to_fit();
//
//    for (const auto & it : variance.delMap) {
//        //std::string str = "deletion: " + std::to_string(it.first) + " -> " + it.second;
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "deletion: " << it.first << " -> " << it.second << "\n";
//    }
//
//    //get<2>(*vTuple).clear();
//    //get<2>(*vTuple).shrink_to_fit();
//
//
//    int tIdx = -1;
//    int qIdx = -1;
//    if (mOptions->mLocSnps.mLocSnpOptions.modeCode == EDLIB_MODE_HW) {
//        tIdx = *(result.endLocations);
//        for (int i = 0; i < result.alignmentLength; i++) {
//            if (result.alignment[i] != EDLIB_EDOP_INSERT) {
//                tIdx--;
//            }
//        }
//    }
//
//    for (int start = 0; start < result.alignmentLength; start += 150) {
//        // target
//        //printf("T: ");
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "T: ";
//        int startTIdx = -1;
//        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
//            if (result.alignment[j] == EDLIB_EDOP_INSERT) {
//                //printf("-");
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << "-";
//            } else {
//                //printf("%c", target[++tIdx]);
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << tData[++tIdx];
//            }
//            if (j == start) {
//                startTIdx = tIdx;
//            }
//        }
//        //printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);
//        if (mOptions->mLocVars.locVarOptions.printRes) {
//            ss << " (" << max(startTIdx, 0) << " - " << tIdx << ")\n";
//            ss << "   ";
//        }
//        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
//            //printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
//            if (mOptions->mLocVars.locVarOptions.printRes) ss << (result.alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
//        }
//        //printf("\n");
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "\n";
//
//        // query
//        //printf("Q: ");
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << "Q: ";
//        int startQIdx = qIdx;
//        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
//            if (result.alignment[j] == EDLIB_EDOP_DELETE) {
//                //printf("-");
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << "-";
//            } else {
//                //printf("%c", query[++qIdx]);
//                if (mOptions->mLocVars.locVarOptions.printRes) ss << qData[++qIdx];
//            }
//            if (j == start) {
//                startQIdx = qIdx;
//            }
//        }
//        //printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
//        if (mOptions->mLocVars.locVarOptions.printRes) ss << " (" << max(startQIdx, 0) << " - " << qIdx << ")\n";
//    }
//}

std::tuple<int, int, bool> SnpScanner::doPrimerAlignment(const char* & qData, int qLength, const std::string & qName,
        const char* & tData, int tLength, const std::string & tName, bool printAlignment) {

    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {

        std::set<int> snpsSet;
        std::set<int> indelSet;
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsSet.insert(i);
            } else if (cur == EDLIB_EDOP_INSERT) {
                indelSet.insert(i);
            } else if (cur == EDLIB_EDOP_DELETE) {
                indelSet.insert(i);
            }
        }

        int endPos = *(result.endLocations) + 1;
        edlibFreeAlignResult(result);
        if (indelSet.empty() && snpsSet.size() <= mOptions->mLocVars.locVarOptions.maxMismatchesPSeq) {
            return std::make_tuple(snpsSet.size(), endPos, true);
        } else {
            return std::make_tuple(0, 0, false);
        }
    } else {
        edlibFreeAlignResult(result);
        return std::make_tuple(0, 0, false);
    }
}

//void SnpScanner::merge(Options * & mOptions, std::vector<std::map<std::string, std::map<std::string, uint32>>> & totalSnpSeqMapVec) {
//    //std::map<std::string, std::map<std::string, LocSnp>> allGenotypeSnpMap;
//    if (totalSnpSeqMapVec.empty()) {
//        return;
//    }
//
//    std::map<std::string, std::map < std::string, uint32>> tmpSnpSeqsMap;
//    for (const auto & it : totalSnpSeqMapVec) {
//        for (const auto & it2 : it) {
//            for (const auto & it3 : it2.second) {
//                tmpSnpSeqsMap[it2.first][it3.first] += (mOptions->isPaired() ? (2 * it3.second) : it3.second);
//            }
//        }
//    }
//
//    std::string foutName = mOptions->prefix + "_snps_genotypes.txt";
//    std::ofstream * fout = new std::ofstream();
//    fout->open(foutName.c_str(), std::ofstream::out);
//
//    if (!fout->is_open()) {
//        delete fout;
//        fout = nullptr;
//        error_exit("Can not open output file: " + foutName);
//    }
//    if (mOptions->verbose) loginfo("Starting to write genotype table!");
//    *fout << "#Locus\tPosition\tGenotype\tNumReads\tReadsRatio\tTotalReads\tNewSnp\n";
//
//    std::string foutName2 = mOptions->prefix + "_snps_haplotype.txt";
//    std::ofstream * fout2 = new std::ofstream();
//    fout2->open(foutName2.c_str(), std::ofstream::out);
//
//    if (!fout2->is_open()) {
//        delete fout2;
//        fout2 = nullptr;
//        error_exit("Can not open output file: " + foutName2);
//    }
//    if (mOptions->verbose) loginfo("Starting to write haplotype table!");
//
//    *fout2 << "#Locus\tHaplotype\tNumReads\tReadsRatio\tTotalHaploReads\tTotalReads\tConclusive\tIndel\tMicroHaplotype\n";
//
//
//    std::string foutName3 = mOptions->prefix + "_all_amplicon.txt";
//    std::ofstream * fout3 = new std::ofstream();
//    fout3->open(foutName3.c_str(), std::ofstream::out);
//    *fout3 << "#Locus\tAmplicon\tNumReads\tTotalReads\n";
//
//    std::string foutName4 = mOptions->prefix + "_error_rate.txt";
//    std::ofstream * fout4 = new std::ofstream();
//    fout4->open(foutName4.c_str(), std::ofstream::out);
//    *fout4 << "#Locus\tErrorRate\tTotalEffectiveReads\n";
//
//    for (const auto & it : tmpSnpSeqsMap) {
//        if (mOptions->mLocSnps.refLocMap.find(it.first) == mOptions->mLocSnps.refLocMap.end()) {
//            continue;
//        }
//
//        LocSnp2* locSnpIt = &(mOptions->mLocSnps.refLocMap[it.first]);
//        std::map<int, std::map<char, int>> baseFreqMap;
//        //std::set<int> posSet;
//
//        for (const auto & it2 : it.second) {
//            locSnpIt->totReads += it2.second;
//            if (it2.second > locSnpIt->maxReads) {
//                locSnpIt->maxReads = it2.second;
//            }
//        }
//
//        if (locSnpIt->maxReads < mOptions->mLocSnps.mLocSnpOptions.minSeqs) continue;
//
//        auto twoPeaks = getTop2MaxKeyValueVec(it.second);
//        locSnpIt->totHaploReads = twoPeaks.front().second;
//
//        if (twoPeaks.size() == 1) {//one homo allele
//            locSnpIt->ratioHaplo = 1;
//            locSnpIt->genoStr3 = "homo";
//        } else if (twoPeaks.size() == 2) {//two alleles
//            const char* rchar1 = twoPeaks.front().first.c_str();
//            const char* rchar2 = twoPeaks.back().first.c_str();
//            auto mapPair = doAlignment(mOptions, "read1", rchar1, twoPeaks.front().first.length(),
//                    "read2", rchar2, twoPeaks.back().first.length());
//
//            if (mapPair.first) {
//                if (mapPair.second.size() < 2) {
//                    mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerL;
//                } else {
//                    mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerH;
//                }
//            } else {
//                mOptions->mLocSnps.mLocSnpOptions.hmPer = mOptions->mLocSnps.mLocSnpOptions.hmPerH;
//            }
//
//            locSnpIt->ratioHaplo = double(twoPeaks.front().second) / (twoPeaks.front().second + twoPeaks.back().second);
//
//            if (locSnpIt->ratioHaplo >= mOptions->mLocSnps.mLocSnpOptions.hmPer) {//homo
//                twoPeaks.pop_back();
//                twoPeaks.shrink_to_fit();
//                locSnpIt->genoStr3 = "homo"; //also include if it is heter against the ref, eg, ref: AA, target: CC;
//                locSnpIt->ratioHaplo = 1; // if the ratioHaplo for homo is not 1, it should have the seq erros
//            } else if (abs(locSnpIt->ratioHaplo - 0.5) <= mOptions->mLocSnps.mLocSnpOptions.htJetter) {//heter
//                locSnpIt->genoStr3 = "heter";
//                locSnpIt->totHaploReads += twoPeaks.back().second;
//            } else {//inconclusive;
//                locSnpIt->genoStr3 = "inconclusive";
//                locSnpIt->totHaploReads += twoPeaks.back().second;
//            }
//
//        }
//
//        bool indel = false;
//        std::set<int> posCorr; //for correction if there are >= 2 snps and ratio is between 0.65 - 0.90 and inconclusive one
//        for (const auto & it2 : it.second) {
//            bool go = false;
//            if (locSnpIt->maxReads >= mOptions->mLocSnps.mLocSnpOptions.minReads4Filter) {//50
//                if (it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqs) {//5
//                    go = true;
//                } else {
//                    go = false;
//                }
//            } else {
//                if (locSnpIt->genoStr3 == "homo") {
//                    if (it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqs) {//5
//                        go = true;
//                    } else {
//                        go = false;
//                    }
//                } else {
//                    if (it2.second >= mOptions->mLocSnps.mLocSnpOptions.minSeqsPer * locSnpIt->maxReads) {//10%
//                        go = true;
//                    } else {
//                        go = false;
//                    }
//                }
//            }
//
//            if (!go) continue;
//
//            locSnpIt->totEffectReads += it2.second;
//            *fout3 << it.first << "\t" << it2.first << "\t" << it2.second << "\t" << locSnpIt->totReads << "\n";
//
//            const char* target = locSnpIt->ref.mStr.c_str();
//            int targetLength = locSnpIt->ref.length();
//            const char* readSeq = it2.first.c_str();
//            int readLength = it2.first.length();
//
//            SimSnps tmpSimSnps;
//            tmpSimSnps.numReads = it2.second;
//            tmpSimSnps.snpPosSet = locSnpIt->refSnpPosSet;
//
//            auto mapPair = doAlignment(mOptions, "read", readSeq, readLength, locSnpIt->name, target, targetLength);
//
//            if (mapPair.first) {//if there is no indel, that's true, including no snps
//
//                //for sequence error rate;
//                for (int pos = 0; pos < it2.first.length(); pos++) {
//                    baseFreqMap[pos][it2.first[pos]] += it2.second;
//                }
//
//                if (!mapPair.second.empty()) {
//                    for (auto & it3 : mapPair.second) {
//                        tmpSimSnps.snpPosSet.insert(it3.first);
//                        locSnpIt->snpPosSet.insert(it3.first);
//                    }
//                }
//
//                for (const auto & it3 : tmpSimSnps.snpPosSet) {
//                    tmpSimSnps.snpsStr.append(std::to_string(it3 + locSnpIt->trimPos.first));
//                    tmpSimSnps.snpsStr.push_back('(');
//                    tmpSimSnps.snpsStr.push_back(locSnpIt->ref.mStr[it3]);
//                    tmpSimSnps.snpsStr.push_back('|');
//                    tmpSimSnps.snpsStr.push_back(it2.first[it3]);
//                    tmpSimSnps.snpsStr.push_back(')');
//                }
//
//                if (locSnpIt->genoStr3 == "homo") {
//                    if (it2.first == twoPeaks.front().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "homo";
//                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                        locSnpIt->snpPosSetTrueHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                    } else {
//                        tmpSimSnps.genoStr8 = "seqerr";
//                        tmpSimSnps.isHaplo = false;
//                    }
//
//                } else if (locSnpIt->genoStr3 == "heter") {
//                    if (it2.first == twoPeaks.front().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "heter1";
//                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                        locSnpIt->snpPosSetTrueHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                    } else if (it2.first == twoPeaks.back().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "heter2";
//                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                        locSnpIt->snpPosSetTrueHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                    } else {
//                        tmpSimSnps.genoStr8 = "seqerr";
//                        tmpSimSnps.isHaplo = false;
//                    }
//                } else {
//
//                    if (it2.first == twoPeaks.front().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "inHeter1";
//                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                        posCorr.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                    } else if (it2.first == twoPeaks.back().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "inHeter2";
//                        locSnpIt->snpPosSetHaplo.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                        posCorr.insert(tmpSimSnps.snpPosSet.begin(), tmpSimSnps.snpPosSet.end());
//                    } else {
//                        tmpSimSnps.genoStr8 = "seqerr";
//                        tmpSimSnps.isHaplo = false;
//                    }
//
//                }
//                locSnpIt->genoMap[it2.first] = tmpSimSnps;
//            } else {//for indels
//                if (locSnpIt->genoStr3 == "homo") {
//                    if (it2.first == twoPeaks.front().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "indel1";
//                        indel = indel || true;
//                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
//                    }
//                } else if (locSnpIt->genoStr3 == "heter") {
//                    if (it2.first == twoPeaks.front().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "indel1";
//                        indel = indel || true;
//                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
//                    } else if (it2.first == twoPeaks.back().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "indel2";
//                        indel = indel || true;
//                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
//                    }
//                } else {
//                    if (it2.first == twoPeaks.front().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "indel1";
//                        indel = indel || true;
//                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
//                    } else if (it2.first == twoPeaks.back().first) {
//                        tmpSimSnps.isHaplo = true;
//                        tmpSimSnps.genoStr8 = "indel2";
//                        indel = indel || true;
//                        locSnpIt->genoMap[it2.first] = tmpSimSnps;
//                    }
//                }
//            }
//        }
//
//        if (locSnpIt->genoMap.empty()) continue;
//
//        if (indel) {
//            locSnpIt->isIndel = true;
//        } else {
//            if (locSnpIt->genoStr3 == "inconclusive" && (!posCorr.empty())) {
//                int numSnps = 0;
//                for (const auto & it2 : posCorr) {
//                    if (twoPeaks.front().first[it2] != twoPeaks.back().first[it2]) {
//                        numSnps++;
//                    }
//                }
//
//                if (numSnps > 1) {
//                    locSnpIt->snpPosSetTrueHaplo.insert(posCorr.begin(), posCorr.end());
//                    locSnpIt->genoMap[twoPeaks.front().first].isHaplo = true;
//                    locSnpIt->genoMap[twoPeaks.front().first].genoStr8 = "heter1";
//                    locSnpIt->genoMap[twoPeaks.back().first].isHaplo = true;
//                    locSnpIt->genoMap[twoPeaks.back().first].genoStr8 = "heter2";
//                    locSnpIt->genoStr3 = "heter";
//                }
//            }
//        }
//
//        //std::set<int> totPosSet;//get total snps positions only for these true haplotypes also include the inconclusive ones
//
//        locSnpIt->totPosSet.insert(locSnpIt->snpPosSet.begin(), locSnpIt->snpPosSet.end());
//
//        std::string haploStr1 = "";
//        std::string haploStr2 = "";
//
//        if (!locSnpIt->isIndel) {
//            for (const auto & it2 : locSnpIt->snpPosSetHaplo) {
//                bool isPrint = false;
//                SimSnp tmpSimSnp;
//                if (twoPeaks.size() == 1) {
//                    tmpSimSnp.snp1 = tmpSimSnp.snp2 = twoPeaks.front().first[it2];
//                    tmpSimSnp.reads1 = tmpSimSnp.reads2 = (twoPeaks.front().second / 2);
//                    tmpSimSnp.ratio = 1.0;
//                    tmpSimSnp.genoStr3 = "homo";
//                    if (twoPeaks.front().first[it2] == locSnpIt->ref.mStr[it2]) {
//                        tmpSimSnp.color = "green";
//                    } else {
//                        tmpSimSnp.color = "orange";
//                    }
//                    isPrint = true;
//
//                    haploStr1.push_back(twoPeaks.front().first[it2]);
//                    haploStr2.push_back(twoPeaks.front().first[it2]);
//                } else if (twoPeaks.size() == 2) {
//                    if (twoPeaks.front().first[it2] == twoPeaks.back().first[it2]) {//AA, CC
//                        tmpSimSnp.genoStr3 = "homo";
//                        tmpSimSnp.reads1 = tmpSimSnp.reads2 = ((twoPeaks.front().second + twoPeaks.back().second) / 2);
//                        tmpSimSnp.ratio = 1.0;
//                        tmpSimSnp.snp1 = tmpSimSnp.snp2 = twoPeaks.front().first[it2];
//
//                        if (twoPeaks.front().first[it2] == locSnpIt->ref.mStr[it2]) {//AA
//                            tmpSimSnp.color = "green";
//                        } else {//CC
//                            tmpSimSnp.color = "orange";
//                        }
//                        isPrint = true;
//
//                        haploStr1.push_back(twoPeaks.front().first[it2]);
//                        haploStr2.push_back(twoPeaks.back().first[it2]);
//
//                    } else {//AC; CA; CT; 
//                        if (twoPeaks.back().first[it2] == locSnpIt->ref.mStr[it2]) {//CA;
//                            tmpSimSnp.reads1 = twoPeaks.back().second;
//                            tmpSimSnp.reads2 = twoPeaks.front().second;
//                            tmpSimSnp.ratio = double(twoPeaks.back().second) / (twoPeaks.front().second + twoPeaks.back().second);
//                            tmpSimSnp.snp1 = twoPeaks.back().first[it2];
//                            tmpSimSnp.snp2 = twoPeaks.front().first[it2];
//
//                        } else {//AC; CT;
//                            tmpSimSnp.reads1 = twoPeaks.front().second;
//                            tmpSimSnp.reads2 = twoPeaks.back().second;
//                            tmpSimSnp.ratio = double(twoPeaks.front().second) / (twoPeaks.front().second + twoPeaks.back().second);
//                            tmpSimSnp.snp1 = twoPeaks.front().first[it2];
//                            tmpSimSnp.snp2 = twoPeaks.back().first[it2];
//                        }
//
//                        if (locSnpIt->genoStr3 == "heter") {
//                            tmpSimSnp.color = (locSnpIt->refSnpPosSet.find(it2) == locSnpIt->refSnpPosSet.end()) ? "orange" : "red";
//                            tmpSimSnp.genoStr3 = "heter";
//                            isPrint = true;
//                            haploStr1.push_back(twoPeaks.front().first[it2]);
//                            haploStr2.push_back(twoPeaks.back().first[it2]);
//                        } else if (locSnpIt->genoStr3 == "inconclusive") {
//                            tmpSimSnp.color = "transparent";
//                            tmpSimSnp.genoStr3 = "inconclusive";
//                            isPrint = false;
//                        }
//                    }
//                }
//                locSnpIt->snpsMap[it2] = tmpSimSnp;
//                if (isPrint) *fout << it.first << "\t" << (it2 + locSnpIt->trimPos.first) << "\t" << tmpSimSnp.snp1 << "|" << tmpSimSnp.snp2 << "\t" << tmpSimSnp.reads1 << "|" << tmpSimSnp.reads2 << "\t" << tmpSimSnp.ratio << "\t" << (tmpSimSnp.reads1 + tmpSimSnp.reads2) << "\t" << (locSnpIt->refSnpPosSet.find(it2) == locSnpIt->refSnpPosSet.end() ? "Y" : "N") << "\n";
//            }
//        }
//
//        if (locSnpIt->isIndel) {
//            if (twoPeaks.size() == 1) {
//                if (locSnpIt->genoStr3 == "inconclusive") {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'Y'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'Y'));
//                } else {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'Y'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'Y'));
//                }
//            } else if (twoPeaks.size() == 2) {
//                if (locSnpIt->genoStr3 == "inconclusive") {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'N', 'Y'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'N', 'Y'));
//                } else {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'Y', 'Y'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'Y', 'Y'));
//                }
//            }
//
//        } else {
//
//            if (twoPeaks.size() == 1) {
//                if (locSnpIt->genoStr3 == "inconclusive") {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'N'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'N', 'N'));
//                } else {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'N'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, (twoPeaks.front().second / 2), (double((twoPeaks.front().second / 2)) / locSnpIt->totHaploReads), 'Y', 'N'));
//
//                    for (int pos = 0; pos < twoPeaks.front().first.length(); pos++) {
//                        int posReads = 0;
//                        auto posIt = baseFreqMap[pos].find(twoPeaks.front().first[pos]);
//                        if (posIt != baseFreqMap[pos].end()) {
//                            for (const auto & posIt2 : baseFreqMap[pos]) {
//                                if (posIt2.first != twoPeaks.front().first[pos]) {
//                                    posReads += posIt2.second;
//                                }
//                            }
//                        }
//                        locSnpIt->baseErrorMap[pos] = static_cast<double> (posReads * 100) / static_cast<double> (locSnpIt->totEffectReads);
//                    }
//
//                }
//            } else if (twoPeaks.size() == 2) {
//                if (locSnpIt->genoStr3 == "inconclusive") {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'N', 'N'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'N', 'N'));
//                } else {
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.front().first, haploStr1, twoPeaks.front().second, (double(twoPeaks.front().second) / locSnpIt->totHaploReads), 'Y', 'N'));
//                    locSnpIt->haploVec.push_back(std::make_tuple(twoPeaks.back().first, haploStr2, twoPeaks.back().second, (double(twoPeaks.back().second) / locSnpIt->totHaploReads), 'Y', 'N'));
//
//                    for (int pos = 0; pos < twoPeaks.front().first.length(); pos++) {
//                        int posReads = 0;
//
//                        if (twoPeaks.front().first[pos] == twoPeaks.back().first[pos]) {
//                            auto posIt = baseFreqMap[pos].find(twoPeaks.front().first[pos]);
//                            if (posIt != baseFreqMap[pos].end()) {
//                                for (const auto & posIt2 : baseFreqMap[pos]) {
//                                    if (posIt2.first != twoPeaks.front().first[pos]) {
//                                        posReads += posIt2.second;
//                                    }
//                                }
//                            }
//                        } else {
//                            auto posIt = baseFreqMap[pos].find(twoPeaks.front().first[pos]);
//                            auto posItb = baseFreqMap[pos].find(twoPeaks.back().first[pos]);
//                            if (posIt != baseFreqMap[pos].end() && posItb != baseFreqMap[pos].end()) {
//                                for (const auto & posIt2 : baseFreqMap[pos]) {
//                                    if (posIt2.first != twoPeaks.front().first[pos] && posIt2.first != twoPeaks.back().first[pos]) {
//                                        posReads += posIt2.second;
//                                    }
//                                }
//                            }
//
//                        }
//
//                        locSnpIt->baseErrorMap[pos] = static_cast<double> (posReads * 100) / static_cast<double> (locSnpIt->totEffectReads);
//                    }
//
//                }
//            }
//
//        }
//
//        if (!locSnpIt->baseErrorMap.empty()) {
//            *fout4 << it.first << "\t";
//            for (const auto & posIt : locSnpIt->baseErrorMap) {
//                if (posIt.first == locSnpIt->baseErrorMap.rbegin()->first) {
//                    *fout4 << posIt.second;
//                } else {
//                    *fout4 << posIt.second << ";";
//                }
//
//            }
//            *fout4 << "\t" << locSnpIt->totEffectReads << "\n";
//        }
//
//        for (const auto & it2 : locSnpIt->haploVec) {
//            locSnpIt->genoMap[get<0>(it2)].haploStr = get<1>(it2);
//            *fout2 << it.first << "\t" << get<1>(it2) << "\t" << get<2>(it2) << "\t" << get<3>(it2) << "\t" << locSnpIt->totHaploReads << "\t" << locSnpIt->totReads << "\t" << get<4>(it2) << "\t" << get<5>(it2) << "\t" << get<0>(it2) << "\n";
//        }
//    }
//
//    tmpSnpSeqsMap.clear();
//
//    fout->flush();
//    fout->clear();
//    fout->close();
//    if (fout) {
//        delete fout;
//        fout = nullptr;
//    }
//
//    if (mOptions->verbose) loginfo("Finished writing genotype table!");
//
//    fout2->flush();
//    fout2->clear();
//    fout2->close();
//    if (fout2) {
//        delete fout2;
//        fout2 = nullptr;
//    }
//
//    if (mOptions->verbose) loginfo("Finished writing haplotype table!");
//
//    fout3->flush();
//    fout3->clear();
//    fout3->close();
//    if (fout3 != nullptr) {
//        delete fout3;
//        fout3 = nullptr;
//    }
//    if (mOptions->verbose) loginfo("Finished writing amplicon table!");
//
//    fout4->flush();
//    fout4->clear();
//    fout4->close();
//    if (fout4 != nullptr) {
//        delete fout4;
//        fout4 = nullptr;
//    }
//    if (mOptions->verbose) loginfo("Finished writing error rate table!");
//
//    //return allGenotypeSnpMap;
//}

void SnpScanner::merge2(Options *&mOptions, std::vector<std::map<std::string, std::map<std::string, uint32>>> &totalSnpSeqMapVec) {

    if (totalSnpSeqMapVec.empty()) {
        return;
    }

    std::map<std::string, std::map<std::string, uint32>> tmpSnpSeqsMap; //merged all seqs;
    for (const auto &it : totalSnpSeqMapVec) {
        for (const auto &it2 : it) {
            for (const auto &it3 : it2.second) {
                tmpSnpSeqsMap[it2.first][it3.first] += (mOptions->isPaired() ? (2 * it3.second) : it3.second);
            }
        }
    }

    std::string foutName = mOptions->prefix + "_snps_genotypes.txt";
    std::ofstream *fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);
    if (!fout->is_open()) {
        delete fout;
        fout = nullptr;
        error_exit("Can not open output file: " + foutName);
    }
    if (mOptions->verbose)
        loginfo("Starting to write snps table!");
    *fout << "#Locus\tPosition\tGenotype\tNumReads\tReadsRatio\tTotalReads\tNewSnp\tConclusive\n";

    std::string foutName2 = mOptions->prefix + "_snps_haplotype.txt";
    std::ofstream *fout2 = new std::ofstream();
    fout2->open(foutName2.c_str(), std::ofstream::out);
    if (!fout2->is_open()) {
        delete fout2;
        fout2 = nullptr;
        error_exit("Can not open output file: " + foutName2);
    }
    if (mOptions->verbose)
        loginfo("Starting to write haplotype table!");

    *fout2 << "#Locus\tHaplotype\tNumHaploReads\tHaploReadsRatio\tHaploReadsPer\tTotalReads\tZegosity\tIndel\tMicroHaplotype\n";

    std::string foutName3 = mOptions->prefix + "_all_amplicon.txt";
    std::ofstream *fout3 = new std::ofstream();
    fout3->open(foutName3.c_str(), std::ofstream::out);
    if (!fout3->is_open()) {
        delete fout3;
        fout3 = nullptr;
        error_exit("Can not open output file: " + foutName3);
    }
    if (mOptions->verbose)
        loginfo("Starting to write haplotype table!");
    *fout3 << "#Locus\tNumReads\tTotalReads\tReadRatio\tSNP\tAmplicon\n";

    std::string foutName4 = mOptions->prefix + "_error_rate.txt";
    std::ofstream *fout4 = new std::ofstream();
    fout4->open(foutName4.c_str(), std::ofstream::out);
    if (!fout4->is_open()) {
        delete fout2;
        fout2 = nullptr;
        error_exit("Can not open output file: " + foutName4);
    }
    if (mOptions->verbose)
        loginfo("Starting to write error rate table!");
    *fout4 << "#Locus\tErrorRate\tTotalReads\n";

    for (const auto &it : tmpSnpSeqsMap) {
        if (mOptions->mLocSnps.refLocMap.find(it.first) == mOptions->mLocSnps.refLocMap.end()) {
            continue;
        }
        if (it.second.empty()) continue;
        LocSnp2 *locSnpIt = &(mOptions->mLocSnps.refLocMap[it.first]);

        std::map<std::string, uint32> tmpMap;
        if (locSnpIt->trimPos.first == 0 && locSnpIt->trimPos.second == 0) {
            tmpMap = it.second;
        } else {
            for (const auto it2 : it.second) {
                std::string str = "";
                bool trimmed = false;
                if (locSnpIt->trimPos.first != 0) {
                    if ((locSnpIt->trimPos.first + locSnpIt->ref.length() + locSnpIt->trimPos.second) <= it2.first.length()) {
                        str = it2.first.substr(locSnpIt->trimPos.first);
                        trimmed = true;
                    } else {
                        continue;
                    }
                } else {
                    str = it2.first;
                }
                
                if(locSnpIt->trimPos.second != 0){
                    if(trimmed){
                        if ((locSnpIt->ref.length() + locSnpIt->trimPos.second) <= str.length()) {
                            str = str.substr(0, str.length() - locSnpIt->trimPos.second);
                        } else {
                            continue;
                        }
                    } else {
                        if ((locSnpIt->trimPos.first + locSnpIt->ref.length() + locSnpIt->trimPos.second) <= str.length()) {
                            str = str.substr(0, str.length() - locSnpIt->trimPos.second);
                        } else {
                            continue;
                        }
                    }
                }
                
                tmpMap[str] += it2.second;
            }
        }
        
        if(tmpMap.empty()) continue;
        
        locSnpIt->seqVarVec.reserve(tmpMap.size());
        std::map<int, std::map<char, int>> baseFreqMap;
        for (const auto &it2 : tmpMap) {
            locSnpIt->totReads += it2.second;
            if (it2.second > locSnpIt->maxReads) {
                locSnpIt->maxReads = it2.second;
            }
            
            const char *target = locSnpIt->ref.mStr.c_str();
            int targetLength = locSnpIt->ref.length();
            const char *readSeq = it2.first.c_str();
            int readLength = it2.first.length();
            auto mapPair = doAlignment2(mOptions, "read", readSeq, readLength, locSnpIt->name, target, targetLength);

            SeqVar tmpSeqVar;
            tmpSeqVar.seq = it2.first;
            tmpSeqVar.numReads = it2.second;
            if (mapPair.first) { // no indel
                tmpSeqVar.snpSet = mapPair.second;
                for (int i = 0; i < it2.first.length(); i++) {
                    baseFreqMap[i][it2.first[i]] += it2.second;
                }
            } else { // has indel;
                tmpSeqVar.indel = true;
            }
            locSnpIt->seqVarVec.push_back(tmpSeqVar);
        }
        
        if (locSnpIt->maxReads < mOptions->mLocSnps.mLocSnpOptions.minSeqs) {
            locSnpIt->seqVarVec.clear();
            baseFreqMap.clear();
            continue;
        }
        std::sort(locSnpIt->seqVarVec.begin(), locSnpIt->seqVarVec.end(), [](const SeqVar & L, const SeqVar & R) {return L.numReads > R.numReads;});

        if (locSnpIt->seqVarVec.size() == 1) {
            locSnpIt->genoStr3 = "homo";
            locSnpIt->ratioHaplo = 1;

            locSnpIt->status.second = locSnpIt->seqVarVec.at(0).indel;
            locSnpIt->status.first.first = locSnpIt->seqVarVec.at(0).indel;

            locSnpIt->totHaploReads = locSnpIt->seqVarVec.at(0).numReads;
            if (locSnpIt->status.second) {
                baseFreqMap.clear();
            } else {
                locSnpIt->snpPosSetHaplo = locSnpIt->refSnpPosSet;
                locSnpIt->snpPosSetHaplo.insert(locSnpIt->seqVarVec.at(0).snpSet.begin(), locSnpIt->seqVarVec.at(0).snpSet.end());
                for (int pos = 0; pos < locSnpIt->seqVarVec.at(0).seq.length(); pos++) {
                    baseFreqMap[pos][locSnpIt->seqVarVec.at(0).seq[pos]] = 0;
                }
            }
        } else {
            double ratio = getPer(locSnpIt->seqVarVec.at(0).numReads, (locSnpIt->seqVarVec.at(0).numReads + locSnpIt->seqVarVec.at(1).numReads), false);
            locSnpIt->ratioHaplo = ratio;
            const char *rchar1 = locSnpIt->seqVarVec.at(0).seq.c_str();
            const char *rchar2 = locSnpIt->seqVarVec.at(1).seq.c_str();
            auto mapPair = doAlignment2(mOptions, "read1", rchar1, locSnpIt->seqVarVec.at(0).seq.length(), "read2", rchar2, locSnpIt->seqVarVec.at(1).seq.length());

            if (mapPair.first) {
                if (mapPair.second.size() == 1) {//one snp;
                    if (ratio >= mOptions->mLocSnps.mLocSnpOptions.hmPerL) {
                        locSnpIt->genoStr3 = "homo"; // also include if it is heter against the ref, eg, ref: AA, target: CC;
                        //locSnpIt->ratioHaplo = 1; // if the ratioHaplo for homo is not 1, it should have the seq erros
                        locSnpIt->status.second = locSnpIt->seqVarVec.at(0).indel;
                    } else if (abs(ratio - 0.5) <= mOptions->mLocSnps.mLocSnpOptions.htJetter) {
                        if (locSnpIt->seqVarVec.at(1).numReads < mOptions->mLocSnps.mLocSnpOptions.minSeqs) {
                            locSnpIt->genoStr3 = "inconclusive";
                        } else {
                            locSnpIt->genoStr3 = "heter";
                        }
                        locSnpIt->status.second = !mapPair.first; // false;
                        locSnpIt->status.first.first = locSnpIt->seqVarVec.at(0).indel;
                        locSnpIt->status.first.second = locSnpIt->seqVarVec.at(1).indel;
                    } else {
                        locSnpIt->genoStr3 = "inconclusive";
                        locSnpIt->status.second = !mapPair.first;
                        locSnpIt->status.first.first = locSnpIt->seqVarVec.at(0).indel;
                        locSnpIt->status.first.second = locSnpIt->seqVarVec.at(1).indel;
                    }
                } else if (mapPair.second.size() > 1) { // > 1 snp;
                    if (ratio >= mOptions->mLocSnps.mLocSnpOptions.hmPerH) {
                        locSnpIt->genoStr3 = "homo"; // also include if it is heter against the ref, eg, ref: AA, target: CC;
                        //locSnpIt->ratioHaplo = 1; // if the ratioHaplo for homo is not 1, it should have the seq erros
                        locSnpIt->status.second = locSnpIt->seqVarVec.at(0).indel;
                    } else {
                        if (locSnpIt->seqVarVec.at(1).numReads < mOptions->mLocSnps.mLocSnpOptions.minSeqs) {
                            locSnpIt->genoStr3 = "inconclusive";
                        } else {
                            locSnpIt->genoStr3 = "heter";
                        }
                        locSnpIt->status.second = !mapPair.first;
                        locSnpIt->status.first.first = locSnpIt->seqVarVec.at(0).indel;
                        locSnpIt->status.first.second = locSnpIt->seqVarVec.at(1).indel;
                    }
                }
            } else {
                if (ratio >= mOptions->mLocSnps.mLocSnpOptions.hmPerH) {// 5 / 5 + x >= 0.9
                    locSnpIt->genoStr3 = "homo"; // also include if it is heter against the ref, eg, ref: AA, target: CC;
                    //locSnpIt->ratioHaplo = 1; // if the ratioHaplo for homo is not 1, it should have the seq erros
                    locSnpIt->status.second = locSnpIt->seqVarVec.at(0).indel;
                } else {
                    if (locSnpIt->seqVarVec.at(1).numReads < mOptions->mLocSnps.mLocSnpOptions.minSeqs) {
                        locSnpIt->genoStr3 = "inconclusive";
                    } else {
                        locSnpIt->genoStr3 = "heter";
                    }
                    locSnpIt->status.second = !mapPair.first;
                    locSnpIt->status.first.first = locSnpIt->seqVarVec.at(0).indel;
                    locSnpIt->status.first.second = locSnpIt->seqVarVec.at(1).indel;
                }
            }
            
            if (locSnpIt->genoStr3 == "homo") {
                locSnpIt->totHaploReads = locSnpIt->seqVarVec.at(0).numReads;
                if (locSnpIt->status.first.first) {
                    baseFreqMap.clear();
                } else {
                    locSnpIt->snpPosSetHaplo = locSnpIt->refSnpPosSet;
                    locSnpIt->snpPosSetHaplo.insert(locSnpIt->seqVarVec.at(0).snpSet.begin(), locSnpIt->seqVarVec.at(0).snpSet.end());
                    for (int pos = 0; pos < locSnpIt->seqVarVec.at(0).seq.length(); pos++) {
                        baseFreqMap[pos][locSnpIt->seqVarVec.at(0).seq[pos]] = 0;
                    }
                }
            } else {
                locSnpIt->totHaploReads = locSnpIt->seqVarVec.at(0).numReads + locSnpIt->seqVarVec.at(1).numReads;
                if (locSnpIt->status.second) {
                    if (locSnpIt->status.first.first && locSnpIt->status.first.second) {
                        //h1 Y, h2 Y;
                        baseFreqMap.clear();
                    } else if (locSnpIt->status.first.first && !locSnpIt->status.first.second) {// h1 Y, h2, N
                        baseFreqMap.clear();
                    } else if (!locSnpIt->status.first.first && locSnpIt->status.first.second) {
                        baseFreqMap.clear();
                    } else {

                    }
                } else {
                    if (locSnpIt->status.first.first && locSnpIt->status.first.second) {
                            baseFreqMap.clear();
                    } else if (locSnpIt->status.first.first && !locSnpIt->status.first.second) {// h1 Y, h2, N
                        baseFreqMap.clear();
                    } else if (!locSnpIt->status.first.first && locSnpIt->status.first.second) {
                        baseFreqMap.clear();
                    } else {
                        locSnpIt->snpPosSet.insert(mapPair.second.begin(), mapPair.second.end());
                        locSnpIt->snpPosSetHaplo = locSnpIt->refSnpPosSet;
                        locSnpIt->snpPosSetHaplo.insert(locSnpIt->seqVarVec.at(0).snpSet.begin(), locSnpIt->seqVarVec.at(0).snpSet.end());
                        locSnpIt->snpPosSetHaplo.insert(locSnpIt->seqVarVec.at(1).snpSet.begin(), locSnpIt->seqVarVec.at(1).snpSet.end());
                        locSnpIt->snpPosSetHaplo.insert(mapPair.second.begin(), mapPair.second.end());
                        for (int pos = 0; pos < locSnpIt->seqVarVec.at(0).seq.length(); pos++) {
                            baseFreqMap[pos][locSnpIt->seqVarVec.at(0).seq[pos]] = 0;
                            baseFreqMap[pos][locSnpIt->seqVarVec.at(1).seq[pos]] = 0;
                        }
                    }
                }
            }
        }

        if (locSnpIt->genoStr3 == "homo") {
            if (!locSnpIt->status.first.first) {
                for (const auto & its : locSnpIt->snpPosSetHaplo) {
                    if (locSnpIt->ref.mStr[its] == locSnpIt->seqVarVec.at(0).seq[its]) {
                        locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->ref.mStr[its], locSnpIt->ref.mStr[its], 'g');
                    } else {
                        if (locSnpIt->refSnpPosSet.find(its) == locSnpIt->refSnpPosSet.end()) {
                            locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->ref.mStr[its], locSnpIt->ref.mStr[its], 'o');
                        } else {
                            locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->ref.mStr[its], locSnpIt->ref.mStr[its], 'r');
                        }
                    }
                }
            }
        } else {
            if (!locSnpIt->status.first.first && !locSnpIt->status.first.second && !locSnpIt->status.second) {
                for (const auto & its : locSnpIt->snpPosSetHaplo) {
                    //3 comparisions;
                    if(locSnpIt->seqVarVec.at(0).seq[its] == locSnpIt->seqVarVec.at(1).seq[its]){
                        if(locSnpIt->ref.mStr[its] == locSnpIt->seqVarVec.at(0).seq[its]){
                            locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->ref.mStr[its], locSnpIt->ref.mStr[its], 'g');
                        } else {
                            if (locSnpIt->refSnpPosSet.find(its) == locSnpIt->refSnpPosSet.end()) {
                                locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->seqVarVec.at(0).seq[its], locSnpIt->seqVarVec.at(1).seq[its], 'o');
                            } else {
                                locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->seqVarVec.at(0).seq[its], locSnpIt->seqVarVec.at(1).seq[its], 'r');
                            }
                        }
                    } else {
                        if (locSnpIt->refSnpPosSet.find(its) == locSnpIt->refSnpPosSet.end()) {
                            locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->seqVarVec.at(0).seq[its], locSnpIt->seqVarVec.at(1).seq[its], 'o');
                        } else {
                            locSnpIt->ssnpsMap[its] = SSimSnp(locSnpIt->seqVarVec.at(0).seq[its], locSnpIt->seqVarVec.at(1).seq[its], 'r');
                        }
                    }
                }
            }
        }

        size_t rmPos = 0;
        for (size_t i = 0; i < locSnpIt->seqVarVec.size(); i++) {
            rmPos = i;
            if (locSnpIt->seqVarVec.at(i).indel) {
                continue;
            }

            bool go = false;
            if (i == 0) {
                go = true;
            } else if (i == 1) {
                go = true;
            } else {
                if (locSnpIt->maxReads >= mOptions->mLocSnps.mLocSnpOptions.minReads4Filter) {
                    if (locSnpIt->seqVarVec.at(i).numReads >= std::max(mOptions->mLocSnps.mLocSnpOptions.minSeqs, (uint32)(mOptions->mLocSnps.mLocSnpOptions.minSeqsPer * locSnpIt->maxReads))) {
                        go = true;
                    }
                } else {
                    if (locSnpIt->genoStr3 == "homo") {
                        if (locSnpIt->seqVarVec.at(i).numReads >= std::max(mOptions->mLocSnps.mLocSnpOptions.minSeqs, (uint32)(mOptions->mLocSnps.mLocSnpOptions.minSeqsPer * locSnpIt->maxReads))) {
                            go = true;
                        }
                    } else {
                        go = true;
                    }
                }
            }

            if (go) {
                locSnpIt->snpPosSet.insert(locSnpIt->seqVarVec.at(i).snpSet.begin(), locSnpIt->seqVarVec.at(i).snpSet.end());
            } else {
                break;
            }
        }
        
        locSnpIt->seqVarVec.erase(locSnpIt->seqVarVec.begin() + (rmPos++), locSnpIt->seqVarVec.end());
        locSnpIt->seqVarVec.shrink_to_fit();

        if(!locSnpIt->ssnpsMap.empty()){
            for(const auto & itss : locSnpIt->ssnpsMap){
                *fout << locSnpIt->name << "\t" << (itss.first + locSnpIt->trimPos.first) << "\t" << itss.second.snp1 << "|" << itss.second.snp2 << "\t";
                if (locSnpIt->genoStr3 == "homo") {
                    *fout << locSnpIt->getHaploReads() << "|" << locSnpIt->getHaploReads() << "\t";
                } else {
                    *fout << locSnpIt->getHaploReads() << "|" << locSnpIt->getHaploReads(true) << "\t";
                }
                *fout << locSnpIt->ratioHaplo << "\t" << locSnpIt->totHaploReads << "\t" << 
                        (locSnpIt->refSnpPosSet.find(itss.first) == locSnpIt->refSnpPosSet.end() ? "Y" : "N") << "\t" << (locSnpIt->genoStr3 == "inconclusive" ? "N" : "Y") << "\n";
            }
        }
        
        //for *fout2;
        //*fout2 << "#Locus\tHaplotype\tNumHaploReads\tHaploReadsRatio\tHaploReadsPer\tTotalReads\tZegosity\tIndel\tMicroHaplotype\n";
        if(locSnpIt->genoStr3 == "homo"){
            *fout2 << locSnpIt->name << "\t" <<locSnpIt->getHaploStr() << "\t" << locSnpIt->getHaploReads() << "\t" << locSnpIt->ratioHaplo << "\t" << 
                    locSnpIt->getHaploReadsPer() << "\t" << locSnpIt->totReads << "\t" << locSnpIt->genoStr3 << "\t" << (locSnpIt->status.first.first ? "Y" : "N") << "\t" << locSnpIt->seqVarVec.at(0).seq << "\n";
            
            *fout2 << locSnpIt->name << "\t" <<locSnpIt->getHaploStr() << "\t" << locSnpIt->getHaploReads() << "\t" << locSnpIt->ratioHaplo << "\t" << 
                    locSnpIt->getHaploReadsPer() << "\t" << locSnpIt->totReads << "\t" << locSnpIt->genoStr3 << "\t" << (locSnpIt->status.first.first ? "Y" : "N") << "\t" << locSnpIt->seqVarVec.at(0).seq << "\n";
        } else {
            *fout2 << locSnpIt->name << "\t" << locSnpIt->getHaploStr() << "\t" << locSnpIt->getHaploReads() << "\t" << locSnpIt->ratioHaplo << "\t" << 
                    locSnpIt->getHaploReadsPer() << "\t" << locSnpIt->totReads << "\t" << locSnpIt->genoStr3 << "\t" << (locSnpIt->status.first.first ? "Y" : "N") << "\t" << locSnpIt->seqVarVec.at(0).seq << "\n";
            
            *fout2 << locSnpIt->name << "\t" << locSnpIt->getHaploStr(true) << "\t" << locSnpIt->getHaploReads(true) << "\t" << locSnpIt->ratioHaplo << "\t" << 
                    locSnpIt->getHaploReadsPer(true) << "\t" << locSnpIt->totReads << "\t" << locSnpIt->genoStr3 << "\t" << (locSnpIt->status.first.second ? "Y" : "N") << "\t" << locSnpIt->seqVarVec.at(1).seq << "\n";
        }

        for (int i = 0; i < locSnpIt->seqVarVec.size(); i++) {
            *fout3 << locSnpIt->name << "\t" << locSnpIt->seqVarVec.at(i).numReads << "\t" << locSnpIt->totReads << "\t" << 
                    getPer(locSnpIt->seqVarVec.at(i).numReads, locSnpIt->totReads) << "\t" << locSnpIt->getSnpStr(i) << "\t" << locSnpIt->seqVarVec.at(i).seq << "\n";
        }

        if (!baseFreqMap.empty()) {
            *fout4 << locSnpIt->name << "\t";
            for (const auto & itb : baseFreqMap) {
                int tot = 0;
                for (const auto & itb2 : itb.second) {
                    tot += itb2.second;
                }
                double per = getPer(tot, locSnpIt->totReads);
                locSnpIt->baseErrorMap[itb.first] = per;

                if (&itb == &(*baseFreqMap.rbegin())) {
                    *fout4 << per << "\t";
                } else {
                    *fout4 << per << ";";
                }
            }
            *fout4 << locSnpIt->totReads << "\n";
        }
    }

    fout->flush();
    fout->clear();
    fout->close();
    if (fout) {
        delete fout;
        fout = nullptr;
    }

    if (mOptions->verbose)
        loginfo("Finished writing genotype table!");

    fout2->flush();
    fout2->clear();
    fout2->close();
    if (fout2) {
        delete fout2;
        fout2 = nullptr;
    }

    if (mOptions->verbose)
        loginfo("Finished writing haplotype table!");

    fout3->flush();
    fout3->clear();
    fout3->close();
    if (fout3 != nullptr) {
        delete fout3;
        fout3 = nullptr;
    }
    if (mOptions->verbose)
        loginfo("Finished writing amplicon table!");

    fout4->flush();
    fout4->clear();
    fout4->close();
    if (fout4 != nullptr) {
        delete fout4;
        fout4 = nullptr;
    }
    if (mOptions->verbose)
        loginfo("Finished writing error rate table!");
}