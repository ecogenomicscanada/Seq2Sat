#include <valarray>

#include "ssrscanner.h"
#include "options.h"

SsrScanner::SsrScanner(Options* opt) {
    mOptions = opt;
    target = NULL;
    targetLength = 0;
    readSeq = NULL;
    readLength = 0;
    readName = "";
    fpData = NULL;
    fpLength = 0;
    rpData = NULL;
    rpLength = 0;
    ss.str("");
    ffEndPos = 0;
    rfStartPos = 0;
    readWithoutFF.clear();
    readWithoutFFRF.clear();
    locMra.clear();
    tmpAllGenotypeMap.clear();
    //sortedAllGenotypeMap.clear();
    enhancer = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    returnedlocus.clear();
    locVarIt = nullptr;
    checkLoci = true;
    minReadLen = 0;
    mismachesFF = 0;
    mismachesRF = 0;
}

SsrScanner::~SsrScanner() {

}

bool SsrScanner::scanVar(Read* & r1, Read* & r2) {
    return true;
}

void SsrScanner::resetData() {
    readSeq = NULL;
    readLength = 0;
    target = NULL;
    targetLength = 0;
}

std::string SsrScanner::trim2ends(Read* & tmpRead, int startPos, int endPos) {
    std::string rStr = tmpRead->mSeq.mStr;
    return rStr.substr(startPos + 1, endPos - startPos - 1);

}

int SsrScanner::doSubStrAlignment(const char* & pData, int pLength, const char* & rData, int rLength,
        int & numMismaches, int & numDeletions, int & numInsertions, bool rev) {

    int trimPos = 0;

    EdlibAlignResult result = edlibAlign(pData, pLength, rData, rLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            EDLIB_MODE_SHW,
            EDLIB_TASK_PATH,
            NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {

        if (rev) {

            if (*result.endLocations - *result.startLocations >= pLength - mOptions->mEdOptions.maxDeletionRPrimer) {

                for (int i = 0; i < result.alignmentLength; i++) {
                    auto j = static_cast<unsigned> (result.alignment[i]);
                    if (j == 3) {
                        numMismaches++;
                    } else if (j == 2) {
                        numDeletions++;
                    } else if (j == 1) {
                        numInsertions++;
                    }
                }

                if (numMismaches <= mOptions->mEdOptions.maxMismatchesRPrimer &&
                        numDeletions <= mOptions->mEdOptions.maxDeletionRPrimer &&
                        numInsertions <= mOptions->mEdOptions.maxInsertionRPrimer) {

                    trimPos = rLength - (*result.endLocations + 1);
                    if (mOptions->mEdOptions.printRes) {
                        ss << "rev_RPrimer: " << pData << " (" << pLength << " residues) against rev_read " << readName << ": (" << rLength << " residues)" << ": score = " << result.editDistance << "\n";
                        //cCout();
                        printAlignment(pData, rData, result.alignment, result.alignmentLength, *(result.endLocations), mOptions->mEdOptions.modeCode);
                    }
                }

            }

        } else {
            if (*result.endLocations - *result.startLocations <= pLength + mOptions->mEdOptions.maxDeletionRPrimer &&
                    rLength - mOptions->mEdOptions.maxDeletionRPrimer <= *result.endLocations &&
                    *result.endLocations < rLength) {

                for (int i = 0; i < result.alignmentLength; i++) {
                    auto j = static_cast<unsigned> (result.alignment[i]);
                    if (j == 3) {
                        numMismaches++;
                    } else if (j == 2) {
                        numDeletions++;
                    } else if (j == 1) {
                        numInsertions++;
                    }
                }

                if (numMismaches <= mOptions->mEdOptions.maxMismatchesRPrimer &&
                        numDeletions <= mOptions->mEdOptions.maxDeletionRPrimer &&
                        numInsertions <= mOptions->mEdOptions.maxInsertionRPrimer) {
                    trimPos = *result.startLocations;
                }

                if (mOptions->mEdOptions.printRes) {
                    ss << "RPrimer: " << pData << " (" << pLength << " residues) against read " << readName << ": (" << rLength << " residues)" << ": score = " << result.editDistance << "\n";
                    //cCout();
                    printAlignment(pData, rData, result.alignment, result.alignmentLength, *(result.endLocations), mOptions->mEdOptions.modeCode);
                }
            }
        }
    }
    edlibFreeAlignResult(result);
    return trimPos;
}

std::pair<std::map<int, std::string>, bool> SsrScanner::doSimpleAlignment(Options * & mOptions, const char* & qData, int qLength, const char* & tData, int tLength) {
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            mOptions->mEdOptions.modeCode,
            mOptions->mEdOptions.alignTask,
            NULL, 0));

    std::map<int, std::string> snpsMap;
    std::set<int> indelSet;
    if (result.status == EDLIB_STATUS_OK) {
        for (int i = 0; i < result.alignmentLength; i++) {
            auto cur = result.alignment[i];
            if (cur == EDLIB_EDOP_MATCH) {

            } else if (cur == EDLIB_EDOP_MISMATCH) {
                snpsMap[i] = tData[i];
            } else if (cur == EDLIB_EDOP_INSERT) {
                indelSet.insert(i);
            } else if (cur == EDLIB_EDOP_DELETE) {
                indelSet.insert(i);
            }
        }

    }
    edlibFreeAlignResult(result);
    if (indelSet.empty()) {
        return std::make_pair(snpsMap, true);
    } else {
        if(qLength == tLength){
            return std::make_pair(snpsMap, true);
        } else {
            return std::make_pair(snpsMap, false);
        }
    }
}

void SsrScanner::printAlignment(const char* & query, const char* & target,
        const unsigned char* alignment, const int alignmentLength,
        const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        ss << "T: ";
        int startTIdx = -1;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                ss << "-";
            else
                ss << target[++tIdx];
            if (j == start)
                startTIdx = tIdx;
        }
        ss << " (" << max(startTIdx, 0) << " - " << tIdx << ")\n";

        // match / mismatch
        ss << "   ";
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            ss << (alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        ss << "\n";

        // query
        ss << "Q: ";
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                ss << "-";
            else
                ss << query[++qIdx];
            if (j == start)
                startQIdx = qIdx;
        }
        ss << " (" << max(startQIdx, 0) << " - " << qIdx << ")\n\n";
    }
}

int SsrScanner::doAlignment(const char* & qData, int qLength, const std::string & qName,
        const char* & tData, int tLength, const std::string & tName, Variance & variance, bool printAlignment) {
    int endPos = 0;
    EdlibAlignResult result = edlibAlign(qData, qLength, tData, tLength,
            edlibNewAlignConfig(mOptions->mLocVars.locVarOptions.maxScorePrimer,
            mOptions->mEdOptions.modeCode,
            mOptions->mEdOptions.alignTask,
            NULL, 0));

    if (result.status == EDLIB_STATUS_OK) {
        endPos = *(result.endLocations);
        doScanVariance(result, variance, qData, tData, *(result.endLocations));
        if (mOptions->mEdOptions.printRes && printAlignment) {
            printVariance(result, variance, qData, qName, tData, tName, *(result.endLocations));
        }
    }
    edlibFreeAlignResult(result);
    return endPos;
}

void SsrScanner::doScanVariance(EdlibAlignResult & result, Variance & variance,
        const char* & qData, const char* & tData, const int position) {

    int ti = -1, qi = -1;
    std::string inStr, deStr, snpStr;

    if (result.alignmentLength < 2) {
        if (mOptions->mEdOptions.printRes) ss << "You must have at least 2 bp!\n";
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

void SsrScanner::printVariance(EdlibAlignResult & result, Variance & variance,
        const char* & qData, const std::string & qName, const char* & tData, const std::string & tName, const int position) {

    if (mOptions->mEdOptions.printRes) {
        ss << "End location: " << *result.endLocations << " -> alignmentLength: " << result.alignmentLength << " with alignment distance: " << result.editDistance << "\n";
        ss << "Ref: " << tName << " -> Query: " << qName << "\n";
        ss << "Repeat seq: " << locVarIt->mra.mStr << " | " << locVarIt->mra.reverseComplement().mStr << "\n";
    }
    for (int i = 0; i < result.alignmentLength; i++) {
        if (mOptions->mEdOptions.printRes) ss << static_cast<unsigned> (result.alignment[i]);
    }
    if (mOptions->mEdOptions.printRes) ss << "\n";


    for (const auto & it : variance.subMap) {
        //std::string str = "snp: " + std::to_string(it.first) + " -> " + it.second;
        if (mOptions->mEdOptions.printRes) ss << "snp: " << it.first << " -> " << it.second << "\n";
    }

    //get<0>(*vTuple).clear();
    //get<0>(*vTuple).shrink_to_fit();

    for (const auto & it : variance.insMap) {
        //std::string str = "insertion: " + std::to_string(it.first) + " -> " + it.second;
        if (mOptions->mEdOptions.printRes) ss << "insertion: " << it.first << " -> " << it.second << "\n";
    }

    //get<1>(*vTuple).clear();
    //get<1>(*vTuple).shrink_to_fit();

    for (const auto & it : variance.delMap) {
        //std::string str = "deletion: " + std::to_string(it.first) + " -> " + it.second;
        if (mOptions->mEdOptions.printRes) ss << "deletion: " << it.first << " -> " << it.second << "\n";
    }

    //get<2>(*vTuple).clear();
    //get<2>(*vTuple).shrink_to_fit();


    int tIdx = -1;
    int qIdx = -1;
    if (mOptions->mEdOptions.modeCode == EDLIB_MODE_HW) {
        tIdx = *(result.endLocations);
        for (int i = 0; i < result.alignmentLength; i++) {
            if (result.alignment[i] != EDLIB_EDOP_INSERT) {
                tIdx--;
            }
        }
    }

    for (int start = 0; start < result.alignmentLength; start += 150) {
        // target
        //printf("T: ");
        if (mOptions->mEdOptions.printRes) ss << "T: ";
        int startTIdx = -1;
        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
            if (result.alignment[j] == EDLIB_EDOP_INSERT) {
                //printf("-");
                if (mOptions->mEdOptions.printRes) ss << "-";
            } else {
                //printf("%c", target[++tIdx]);
                if (mOptions->mEdOptions.printRes) ss << tData[++tIdx];
            }
            if (j == start) {
                startTIdx = tIdx;
            }
        }
        //printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);
        if (mOptions->mEdOptions.printRes) {
            ss << " (" << max(startTIdx, 0) << " - " << tIdx << ")\n";
            ss << "   ";
        }
        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
            //printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
            if (mOptions->mEdOptions.printRes) ss << (result.alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        //printf("\n");
        if (mOptions->mEdOptions.printRes) ss << "\n";

        // query
        //printf("Q: ");
        if (mOptions->mEdOptions.printRes) ss << "Q: ";
        int startQIdx = qIdx;
        for (int j = start; j < start + 150 && j < result.alignmentLength; j++) {
            if (result.alignment[j] == EDLIB_EDOP_DELETE) {
                //printf("-");
                if (mOptions->mEdOptions.printRes) ss << "-";
            } else {
                //printf("%c", query[++qIdx]);
                if (mOptions->mEdOptions.printRes) ss << qData[++qIdx];
            }
            if (j == start) {
                startQIdx = qIdx;
            }
        }
        //printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
        if (mOptions->mEdOptions.printRes) ss << " (" << max(startQIdx, 0) << " - " << qIdx << ")\n";
    }
}

std::pair<size_t, int> SsrScanner::analyzeMRA(std::string rawStr, //rawStr must be passed by value not reference
        const std::string & ssr, std::size_t & ffTrimPos, std::size_t & rfTrimPos) {

    if (ffTrimPos + static_cast<std::size_t> (mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length())) > rfTrimPos) {
        ffTrimPos = 0;
        rfTrimPos = rawStr.length() - 1;
    }

    if (rfTrimPos != 0) {
        rawStr = rawStr.substr(0, rfTrimPos);
    }
    rawStr = rawStr.substr(ffTrimPos);

    std::pair<size_t, int> posSSRlenPair(0, 0);
    std::map<size_t, int> posMap; //pos, n. of ssrs
    size_t startPos = rawStr.find(ssr);
    size_t endPos = 0;
    if (startPos == std::string::npos) {
        return posSSRlenPair;
    }
    while (startPos != std::string::npos) {
        endPos = rawStr.find(ssr, startPos + ssr.size());
        if (endPos == std::string::npos) {
            if (posMap.empty()) {
                posMap[startPos] = 1;
            } else {
                auto it = posMap.rbegin();
                if (it->first + it->second * ssr.length() == startPos) {
                    it->second++;
                } else {
                    posMap[startPos] = 1;
                }
            }
            break;
        } else {
            if (endPos - startPos == ssr.size()) {
                if (posMap.empty()) {
                    posMap[startPos]++;
                } else {
                    posMap.rbegin()->second++;
                }

            } else {
                if (posMap.empty()) {
                    posMap[startPos]++;
                } else {
                    posMap.rbegin()->second++;
                }
                if (endPos != std::string::npos) {
                    posMap[endPos] = 0;
                }
            }
        }
        startPos = endPos;
    }

    size_t coreKey = 0;
    int coreCount = 0;
    std::vector<std::pair<size_t, int>> preVec;
    for (const auto & it : posMap) {
        preVec.emplace_back(std::make_pair(it.first, it.second));
        if (it.second > coreCount) {
            coreCount = it.second;
            coreKey = it.first;
        }
    }

    std::map<size_t, int> posMapFil = {
        {coreKey, coreCount}
    };

    int coreIndex = distance(posMap.begin(), posMap.find(coreKey));
    posMap.clear();
    int index = coreIndex;
    while (index > 0) {
        auto cur = preVec.at(index);
        index--;
        auto pre = preVec.at(index);
        if (cur.first - (pre.first + pre.second * ssr.length()) == ssr.length()) {
            posMapFil[pre.first] = pre.second;
        } else {
            break;
        }
    }
    index = coreIndex;

    while (index < preVec.size() - 1) {
        auto cur = preVec.at(index);
        index++;
        auto pos = preVec.at(index);
        if (pos.first - (cur.first + cur.second * ssr.length()) == ssr.length()) {
            posMapFil[pos.first] = pos.second;
        } else {
            break;
        }
    }
    preVec.clear();

    std::vector<std::pair<size_t, int>> posVec;
    for (const auto & it : posMapFil) {
        posVec.emplace_back(std::make_pair(it.first, it.second));
    }

    coreIndex = distance(posMapFil.begin(), posMapFil.find(coreKey));
    std::pair<size_t, int> corePair = posVec.at(coreIndex);
    posSSRlenPair = std::make_pair(corePair.first, corePair.second * ssr.length());
    int posIndex = coreIndex;
    while (posIndex > 0) {
        auto cur = posVec.at(posIndex);
        posIndex--;
        auto pre = posVec.at(posIndex);
        if (cur.first - (pre.first + pre.second * ssr.length()) == ssr.length()) {
            auto testSeq = rawStr.substr(pre.first, cur.first - pre.first);
            auto testSeqPos = locVarIt->ff.mStr.find(testSeq);
            if (testSeqPos == std::string::npos) {
                if (locVarIt->ff.mStr.length() >= testSeq.length()) {
                    unsigned int ed = edit_distance(testSeq, locVarIt->ff.mStr.substr(locVarIt->ff.mStr.length() - testSeq.length(), testSeq.length()));
                    if (ed <= locVarIt->edCutoff) {
                        posSSRlenPair = std::make_pair(pre.first, posSSRlenPair.first - pre.first + posSSRlenPair.second);
                    } else {
                        break;
                    }
                } else {
                    posSSRlenPair = std::make_pair(pre.first, posSSRlenPair.first - pre.first + posSSRlenPair.second);
                }
            } else {
                if (testSeqPos == locVarIt->ff.mStr.length() - testSeq.length()) {
                    break;
                } else {
                    posSSRlenPair = std::make_pair(pre.first, posSSRlenPair.first - pre.first + posSSRlenPair.second);
                }
            }

        } else {
            break;
        }
    }

    posIndex = coreIndex;
    while (posIndex < posVec.size() - 1) {
        auto cur = posVec.at(posIndex);
        posIndex++;
        auto suc = posVec.at(posIndex);
        if (suc.first - (cur.first + cur.second * ssr.length()) == ssr.length()) {
            auto testSeq = rawStr.substr(cur.first + cur.second * ssr.length(),
                    suc.first + suc.second * ssr.length() - (cur.first + cur.second * ssr.length()));
            auto testSeqPos = locVarIt->rf.mStr.find(testSeq);

            if (testSeqPos == std::string::npos) {

                if (locVarIt->rf.mStr.length() >= testSeq.length()) {
                    unsigned int ed = edit_distance(testSeq, locVarIt->rf.mStr.substr(0, testSeq.length()));
                    if (ed <= locVarIt->edCutoff) {
                        posSSRlenPair = std::make_pair(posSSRlenPair.first, suc.first + suc.second * ssr.length() - posSSRlenPair.first);
                    } else {
                        break;
                    }
                } else {
                    posSSRlenPair = std::make_pair(posSSRlenPair.first, suc.first + suc.second * ssr.length() - posSSRlenPair.first);
                }

            } else {
                if (testSeqPos == 0) {
                    break;
                } else {
                    posSSRlenPair = std::make_pair(posSSRlenPair.first, suc.first + suc.second * ssr.length() - posSSRlenPair.first);
                }
            }
        } else {
            break;
        }
    }

    posMapFil.clear();
    posVec.clear();

    if (mOptions->mEdOptions.printRes) {
        std::string msg = "rawStr: " + rawStr;
        ss << msg << "\n";
        //cCout(msg, 'r');
        msg = "mra: " + rawStr.substr(posSSRlenPair.first, posSSRlenPair.second);
        ss << msg << "\n";
        //cCout(msg, 'b');
    }
    //posSSRlenPair = std::make_pair(posSSRlenPair.first );
    if (rfTrimPos == 0) {
        posSSRlenPair = std::make_pair(posSSRlenPair.first + ffTrimPos, posSSRlenPair.second);
    } else {
        posSSRlenPair = std::make_pair(posSSRlenPair.first + ffTrimPos, rfTrimPos - (posSSRlenPair.first + ffTrimPos));
    }

    return posSSRlenPair;
}

void SsrScanner::analyzeMRASub(Genotype & genotype, std::map<int, std::string> & subMap) {
    std::vector<size_t> subMapVec;
    genotype.baseLocVar.mra.mStr = locVarIt->mra.mStr;
    size_t si = 0;

    for (const auto & it : subMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            si = (size_t) (it.first - locVarIt->ff.length());
            genotype.baseLocVar.mra.mStr.replace(si, 1, it.second.substr(2, 1));
            subMapVec.emplace_back(si);
        }
    }

    genotype.baseLocVar.name = locVarIt->name + "_";
    if (!subMapVec.empty()) {
        int count = 0;
        std::string tmpStr = "";
        std::string tmpGenotypeStr = "";
        size_t startPosSi = 0;
        for (int i = 0; i < locVarIt->nRep; i++) {
            startPosSi = (size_t) (i * locVarIt->repuit.length());
            tmpStr = genotype.baseLocVar.mra.mStr.substr(startPosSi, locVarIt->repuit.length());
            if (tmpStr == locVarIt->repuit.mStr) {
                count++;
                tmpGenotypeStr = "(" + locVarIt->repuit.mStr + ")" + std::to_string(count);
            } else {
                genotype.baseLocVar.name.append(tmpGenotypeStr);
                tmpGenotypeStr.clear();
                count = 0;
                genotype.baseLocVar.name.append(tmpStr);
            }
        }
        if (count != 0) {
            genotype.baseLocVar.name.append(tmpGenotypeStr);
        }
    } else {
        genotype.baseLocVar.name.append("(" + locVarIt->repuit.mStr + ")" + std::to_string(locVarIt->nRep));
    }
}

bool SsrScanner::checkMRASub(Variance & fVar) {
    bool shouldAnalysis = false;

    for (const auto & it : fVar.subMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            ss << "sub111: " << it.first << " : " << it.second << "\n";
            shouldAnalysis = true;
        }
    }

    for (const auto & it : fVar.insMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            ss << "insert111: " << it.first << " : " << it.second << "\n";
            shouldAnalysis = false;
            //return shouldAnalysis;
        }
    }

    for (const auto & it : fVar.delMap) {
        if (it.first >= locVarIt->ff.length() && it.first <= locVarIt->ff.length() + locVarIt->mra.length() - 1) {
            shouldAnalysis = false;
            ss << "del111: " << it.first << " : " << it.second << "\n";
            //return shouldAnalysis;
        }
    }

    return shouldAnalysis;
}

std::map<std::string, std::map<std::string, Genotype>> SsrScanner::merge(std::vector<std::map<std::string, std::map<std::string, Genotype>>> & totalGenotypeSsrMapVec) {

    std::map<std::string, std::map < std::string, Genotype>> allGenotypeSsrMap;
    if (totalGenotypeSsrMapVec.empty()) {
        return allGenotypeSsrMap;
    }

    std::set<std::string> locSet;
    for (const auto & it : totalGenotypeSsrMapVec) {//std::vector<std::map<std::string, std::map< std::string, Genotype>>>
        for (const auto & it2 : it) {
            locSet.insert(it2.first);
        }
    }

    for (const auto & it : locSet) {
        std::map<std::string, Genotype> eachGenotypeMap;
        for (const auto & it2 : totalGenotypeSsrMapVec) {
            auto it3 = it2.find(it);
            if (it3 != it2.end()) {
                for (const auto & it4 : it3->second) {
                    auto it5 = eachGenotypeMap.find(it4.first);
                    if (it5 == eachGenotypeMap.end()) {
                        eachGenotypeMap[it4.first] = it4.second;
                    } else {
                        it5->second.numReads += it4.second.numReads;
                    }
                }
            }
        }
        allGenotypeSsrMap[it] = eachGenotypeMap;
    }
    return allGenotypeSsrMap;
}

std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> SsrScanner::report(Options * & mOptions, std::map<std::string, std::map<std::string, Genotype>> &allGenotypeMap) {
  
    std::map<std::string, std::vector<std::pair < std::string, Genotype>>> sortedAllGenotypeMap;
    std::map<std::string, std::vector<std::pair < std::string, Genotype>>> sortedAllGenotypeMraMap; //marker, seq, geno
    std::vector<std::map < std::string, std::vector<std::pair < std::string, Genotype>>>> sortedAllGenotypeMapVec
    {
        sortedAllGenotypeMap, sortedAllGenotypeMraMap
    };

    const char* target;
    int targetLength;
    const char* readSeq;
    int readLength;

    if (mOptions->isPaired()) {
        for (auto & it : allGenotypeMap) {
            for (auto & it2 : it.second) {
                it2.second.numReads *= 2;
            }
        }
    }

    for (auto & it : allGenotypeMap) {//marker, seq, genotypes;
        auto locVarIt = &(mOptions->mLocVars.refLocMap.find(it.first)->second);
        int maxReads = 0;

        std::map<int, int> tmpReadsMap;
        for (const auto & it2 : it.second) {
            tmpReadsMap[it2.first.length()] += it2.second.numReads;
        }

        for (const auto & it2 : tmpReadsMap) {
            if (it2.second > maxReads) {
                maxReads = it2.second;
            }
        }

        tmpReadsMap.clear();

        std::set<int> ffSet;
        std::set<int> rfSet;
        std::vector<std::pair < std::string, Genotype>> tmpVec;
        tmpVec.reserve(it.second.size());
        std::vector<std::pair < std::string, Genotype>> tmpMraVec;
        tmpMraVec.reserve(it.second.size());
        for (auto & it2 : it.second) {
            if (it2.second.numReads >= mOptions->mLocVars.locVarOptions.minSeqs &&
                    (it2.second.numReads * 100 / maxReads) >= mOptions->mLocVars.locVarOptions.minSeqsPer &&
                    it2.second.baseLocVar.mra.mStr.length() >= (it2.second.baseLocVar.repuit.mStr.length() + it2.second.baseLocVar.repuit2.mStr.length()) * mOptions->mLocVars.locVarOptions.minNSSRUnit) {
      
                auto tmpGeno = getGenotype(it2.second.baseLocVar.mra.mStr, locVarIt->repuit.mStr, locVarIt->repuit2.mStr);
                auto maxNRep = getMaxNumRep(tmpGeno);
                if (maxNRep >= mOptions->mLocVars.locVarOptions.minNSSRUnit) {                                   
                    target = locVarIt->ff.mStr.c_str();
                    targetLength = locVarIt->ff.mStr.length();
                    readSeq = it2.second.baseLocVar.ff.mStr.c_str();
                    readLength = it2.second.baseLocVar.ff.mStr.length();                    
                    auto mapSetff = SsrScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                    if (mapSetff.second) {
                        it2.second.baseLocVar.snpsMapff = mapSetff.first;
                        for (auto & it3 : it2.second.baseLocVar.snpsMapff) {
                            ffSet.insert(it3.first);
                        }
                        target = locVarIt->rf.mStr.c_str();
                        targetLength = locVarIt->rf.mStr.length();
                        readSeq = it2.second.baseLocVar.rf.mStr.c_str();
                        readLength = it2.second.baseLocVar.rf.mStr.length();
                        auto mapSetrf = SsrScanner::doSimpleAlignment(mOptions, readSeq, readLength, target, targetLength);
                        if (mapSetrf.second) {
                            it2.second.baseLocVar.snpsMaprf = mapSetrf.first;
                            for (auto & it3 : it2.second.baseLocVar.snpsMaprf) {
                                rfSet.insert(it3.first);
                            }
                            mOptions->mLocVars.refLocMap[it.first] = *locVarIt;
                            it2.second.baseLocVar.mraName = tmpGeno;
                            it2.second.baseLocVar.repuitAll = locVarIt->repuitAll;
                            it2.second.baseLocVar.mraBase = getMraBase(it2.second.baseLocVar.mraName);
                            tmpVec.push_back(std::make_pair(it2.first, it2.second));
                            tmpMraVec.push_back(std::make_pair(it2.first, it2.second));
                        }
                    }
                }
            }
        }
        mOptions->mLocVars.refLocMap[it.first].refSnpsSetffMap[basename(mOptions->prefix)] = ffSet;
        mOptions->mLocVars.refLocMap[it.first].refSnpsSetrfMap[basename(mOptions->prefix)] = rfSet;
        ffSet.clear();
        rfSet.clear();
        std::sort(tmpVec.begin(), tmpVec.end(),
                [](const std::pair<std::string, Genotype> & l, const std::pair<std::string, Genotype> & r) {
                    if (l.second.baseLocVar.effectiveLen != r.second.baseLocVar.effectiveLen) {
                        return l.second.baseLocVar.effectiveLen < r.second.baseLocVar.effectiveLen;
                    }
                    return l.first.length() < r.first.length();
                });
        sortedAllGenotypeMapVec.at(0)[it.first] = tmpVec;
        tmpVec.clear();
        tmpVec.shrink_to_fit();

        std::sort(tmpMraVec.begin(), tmpMraVec.end(),
                [](const std::pair<std::string, Genotype> & l, const std::pair<std::string, Genotype> & r) {
                    if (l.second.baseLocVar.mra.mStr.length() != r.second.baseLocVar.mra.mStr.length()) {
                        return l.second.baseLocVar.mra.mStr.length() < r.second.baseLocVar.mra.mStr.length();
                    }
                    return l.first.length() < r.first.length();
                });

        sortedAllGenotypeMapVec.at(1)[it.first] = tmpMraVec;
        tmpMraVec.clear();
        tmpMraVec.shrink_to_fit();
    }
    allGenotypeMap.clear();

    std::string foutName = mOptions->prefix + "_genotypes_mra.txt";
    std::ofstream* fout = new std::ofstream();
    fout->open(foutName.c_str(), std::ofstream::out);

    if (!fout->is_open()) error_exit("Can not open output file: " + foutName);
    if (mOptions->verbose) loginfo("Starting to write genotype table!");

    *fout << "#Sample\tLocus\tMicrosatellite\tMRABase\tMRAName\tMRASize\tAllele\tNumReads\tPutativeAllele\tFF\tMRA\tRF\tSnpsFF\tSnpsRF\n";
    for (auto & it : sortedAllGenotypeMapVec.at(1)) {//marker, seq, genotype
        auto locVarIt = &(mOptions->mLocVars.refLocMap[it.first]); //marker
        std::map<int, int> genoReadsMap;
        std::map<int, int> genoReadsMap2;
        std::map<int, int> genoReadsMapT;

        std::map<int, int> ffsnpMap;
        std::map<int, int> rfsnpMap;
        std::map<std::string, int> mraBaseMap;
        for (const auto & it2 : it.second) {
            genoReadsMap[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
            mraBaseMap[it2.second.baseLocVar.mraBase]++;
            auto ffset = mOptions->mLocVars.refLocMap[it.first].refSnpsSetffMap[basename(mOptions->prefix)];

            if (!ffset.empty()) {
                if (!it2.second.baseLocVar.snpsMapff.empty()) {
                    for (const auto & snp : it2.second.baseLocVar.snpsMapff) {
                        ffsnpMap[snp.first]++;
                    }
                }
            }

            auto rfset = mOptions->mLocVars.refLocMap[it.first].refSnpsSetrfMap[basename(mOptions->prefix)];

            if (!rfset.empty()) {
                if (!it2.second.baseLocVar.snpsMaprf.empty()) {
                    for (const auto & snp : it2.second.baseLocVar.snpsMaprf) {
                        rfsnpMap[snp.first]++;
                    }
                }
            }
        }

        if (genoReadsMap.size() == 1) {//only one peak;
            genoReadsMapT[genoReadsMap.begin()->first] = genoReadsMap.begin()->second;
        } else if (genoReadsMap.size() > 1) {// > 1 peak;

            std::pair<std::string, double> mraBaseRatio{"", 0};
            if (mraBaseMap.size() > 1) {

                std::set<std::string> keySet;
                for (const auto & mit : mraBaseMap) {
                    keySet.insert(mit.first);
                }

                std::map<std::string, double> keyRatioMap;
                for (const auto & mkey : keySet) {
                    int alle1 = 0;
                    int alle2 = 0;
                    for (const auto & mit : mraBaseMap) {
                        if (mkey == mit.first) {
                            alle1 += mit.second;
                        } else {
                            alle2 += mit.second;
                        }

                        if (std::min(alle1, alle2) != 0) {
                            double ratio = alle1 > alle2 ? ((double) alle1 / alle2) : ((double) alle2 / alle1);
                            if (ratio >= 1 && ratio <= mOptions->mLocVars.locVarOptions.varRatio) {
                                keyRatioMap[mkey] = ratio;
                            }
                        }
                    }
                }
                if (!keyRatioMap.empty()) {
                    mraBaseRatio = getMinKeyValue(keyRatioMap);
                }
            }


            //firstly analyze the genotype based on snps, which is more accurate, i guess

            std::pair<double, int> ratiofp{0, 0}; //ratio, pos
            std::pair<double, int> ratiorp{0, 0};
            if (!ffsnpMap.empty()) {
                double ratio = 0;
                int pos = 0;
                for (const auto & s : ffsnpMap) {
                    auto geno1 = s.second;
                    if (geno1 != it.second.size()) {
                        auto geno2 = it.second.size() - geno1;
                        auto tmpratio = geno1 > geno2 ? ((double) geno1 / geno2) : ((double) geno2 / geno1);
                        if (tmpratio >= 1 && tmpratio <= mOptions->mLocVars.locVarOptions.varRatio) {
                            if (s.first == ffsnpMap.begin()->first) {
                                ratio = tmpratio;
                                pos = s.first;
                            } else {
                                if (ratio > tmpratio) {
                                    ratio = tmpratio;
                                    pos = s.first;
                                }
                            }
                        }
                    }
                }

                if (ratio != 0) ratiofp = std::make_pair(ratio, pos);
            }

            if (!rfsnpMap.empty()) {
                double ratio = 0;
                int pos = 0;
                for (const auto & s : rfsnpMap) {
                    auto geno1 = s.second;
                    if (geno1 != it.second.size()) {
                        auto geno2 = it.second.size() - geno1;
                        auto tmpratio = geno1 > geno2 ? ((double) geno1 / geno2) : ((double) geno2 / geno1);
                        if (tmpratio >= 1 && tmpratio <= mOptions->mLocVars.locVarOptions.varRatio) {
                            if (s.first == rfsnpMap.begin()->first) {
                                ratio = tmpratio;
                                pos = s.first;
                            } else {
                                if (ratio > tmpratio) {
                                    ratio = tmpratio;
                                    pos = s.first;
                                }
                            }
                        }
                    }
                }

                if (ratio != 0) ratiorp = std::make_pair(ratio, pos);
            }

            std::map<int, int> outMap1;
            std::map<int, int> outMap2;

            bool analBase = false;
            bool analFF = false;
            bool analRF = false;
            if (mraBaseRatio.second == 0) {
                if (ratiofp.first != 0 && ratiorp.first != 0) {
                    if (ratiofp.first < ratiorp.first) {
                        analFF = true;
                    } else {
                        analRF = true;
                    }
                } else if (ratiofp.first != 0 && ratiorp.first == 0) {
                    analFF = true;
                } else if (ratiofp.first == 0 && ratiorp.first != 0) {
                    analRF = true;
                } else {

                }
            } else {

                if (ratiofp.first != 0 && ratiorp.first != 0) {
                    double minRatio = std::min(mraBaseRatio.second, std::min(ratiofp.first, ratiorp.first));

                    if (ratiofp.first == minRatio) {
                        analFF = true;
                    } else if (ratiorp.first == minRatio) {
                        analRF = true;
                    } else {
                        analBase = true;
                    }

                } else if (ratiofp.first != 0 && ratiorp.first == 0) {
                    if (ratiofp.first <= mraBaseRatio.second) {
                        analFF = true;
                    } else {
                        analBase = true;
                    }
                } else if (ratiofp.first == 0 && ratiorp.first != 0) {
                    if (ratiorp.first <= mraBaseRatio.second) {
                        analRF = true;
                    } else {
                        analBase = true;
                    }
                } else {
                    analBase = true;
                }

            }

            if (analFF) {
                for (const auto & it2 : it.second) {
                    if (!it2.second.baseLocVar.snpsMapff.empty()) {
                        bool found = false;
                        for (const auto & it3 : it2.second.baseLocVar.snpsMapff) {
                            if (ratiofp.second == it3.first) {
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            outMap1[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                        } else {
                            outMap2[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                        }
                    } else {
                        outMap2[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                    }
                }
            } else if (analRF) {
                for (const auto & it2 : it.second) {
                    if (!it2.second.baseLocVar.snpsMaprf.empty()) {
                        bool found = false;
                        for (const auto & it3 : it2.second.baseLocVar.snpsMaprf) {
                            if (ratiorp.second == it3.first) {
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            outMap1[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                        } else {
                            outMap2[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                        }
                    } else {
                        outMap2[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                    }
                }
            } else if (analBase) {
                for (const auto & it2 : it.second) {
                    if (mraBaseRatio.first == it2.second.baseLocVar.mraBase) {
                        outMap1[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                    } else {
                        outMap2[it2.second.baseLocVar.effectiveLen] += it2.second.numReads;
                    }
                }
            } else {

            }

            bool shouldGo = false;            
            if (!outMap1.empty() && !outMap2.empty()) {

                auto p1 = getMaxKeyValue(outMap1, true);
                outMap1.clear();
                auto p2 = getMaxKeyValue(outMap2, true);
                outMap2.clear();
                if (p1.first == p2.first) {
                    auto ra = (double) std::min(p1.second, p2.second) / std::max(p1.second, p2.second);
                    if (ra >= mOptions->mLocVars.locVarOptions.heterRatio) {
                        genoReadsMapT.insert(p1);
                        genoReadsMapT.insert(p2);
                        shouldGo = false;
                    } else {
                        shouldGo = true;
                    }
                } else {
                    int diff = std::abs(p1.first - p2.first);
                    double ratio = (double) std::min(p1.second, p2.second) / std::max(p1.second, p2.second);
                    if (p1.first > p2.first) {
                        if (p1.second > p2.second) {
                            if (locVarIt->repuitAllLen == 0 || locVarIt->repuitAllLen == locVarIt->repuit.length()) {
                                if (diff <= locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }
                            } else {
                                if (diff <= locVarIt->repuit.length() || diff <= locVarIt->repuit2.length() || diff <= locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length() || diff <= 2 * locVarIt->repuit2.length() || diff <= 2 * locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }
                            }

                        } else {
                            if (locVarIt->repuitAllLen == 0 || locVarIt->repuitAllLen == locVarIt->repuit.length()) {
                                if (diff <= locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 3 * locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }

                            } else {
                                if (diff <= locVarIt->repuit.length() || diff <= locVarIt->repuit2.length() || diff <= locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length() || diff <= 2 * locVarIt->repuit2.length() || diff <= 2 * locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }
                            }

                        }
                    } else {
                        if (p1.second > p2.second) {
                            if (locVarIt->repuitAllLen == 0 || locVarIt->repuitAllLen == locVarIt->repuit.length()) {
                                if (diff <= locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if(diff <= 3 * locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }
                            } else {
                                if (diff <= locVarIt->repuit.length() || diff <= locVarIt->repuit2.length() || diff <= locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length() || diff <= 2 * locVarIt->repuit2.length() || diff <= 2 * locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if(diff <= 3 * locVarIt->repuit.length() || diff <= 3 * locVarIt->repuit2.length() || diff <= 3 * locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 3) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }
                            }

                        } else {
                            if (locVarIt->repuitAllLen == 0 || locVarIt->repuitAllLen == locVarIt->repuit.length()) {
                                if (diff <= locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length()) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }

                            } else {
                                if (diff <= locVarIt->repuit.length() || diff <= locVarIt->repuit2.length() || diff <= locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else if (diff <= 2 * locVarIt->repuit.length() || diff <= 2 * locVarIt->repuit2.length() || diff <= 2 * locVarIt->repuitAllLen) {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                } else {
                                    if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                        genoReadsMapT.insert(p1);
                                        genoReadsMapT.insert(p2);
                                        shouldGo = false;
                                    } else {
                                        shouldGo = true;
                                    }
                                }
                            }

                        }
                    }
                }

            } else {
                shouldGo = true;
            }

            if (shouldGo) {
                //get the break point of genotypes in x
                int pre = 0, cur = 0;
                auto itg = genoReadsMap.begin();
                pre = itg->first;
                std::vector<int> gapVec;
                std::advance(itg, 1);
                while (itg != genoReadsMap.end()) {
                    cur = itg->first;
                    int gap = cur - pre;
                    if (locVarIt->repuit2.length() == 0) {
                        if (gap != locVarIt->repuit.length()) {
                            gapVec.emplace_back(cur);
                        }
                    } else {
                        if (locVarIt->repuit.length() == locVarIt->repuit2.length()) {
                            if (gap != locVarIt->repuit.length()) {
                                gapVec.emplace_back(cur);
                            }
                        } else {
                            auto maxssr = std::max(locVarIt->repuit.length(), locVarIt->repuit2.length());
                            if (gap != locVarIt->repuit.length() && gap != locVarIt->repuit2.length() && gap > maxssr) {
                                gapVec.emplace_back(cur);
                            }
                        }
                    }
                    pre = cur;
                    itg++;
                }

                bool nogap = false;
                if (!gapVec.empty()) {                   
                    std::map<int, int> maxMap;
                    std::map<int, int> outMap;
                    auto range_pre = genoReadsMap.begin();
                    for (const auto & itgv : gapVec) {
                        auto range_cur = genoReadsMap.find(itgv);
                        if (range_cur != genoReadsMap.end()) {
                            outMap.insert(range_pre, range_cur);
                            maxMap.insert(getMaxKeyValue(outMap, true));
                            outMap.clear();
                            range_pre = range_cur;
                        }
                    }

                    if (range_pre != genoReadsMap.end()) {
                        outMap.insert(range_pre, genoReadsMap.end());
                        maxMap.insert(getMaxKeyValue(outMap, true));
                        outMap.clear();
                    }

                    if (maxMap.size() > 2) {
                        genoReadsMapT = getTop2MaxKeyValue(maxMap);
                    } else {
                        genoReadsMapT = maxMap;
                    }
                    maxMap.clear();
                    int diff = genoReadsMapT.rbegin()->first - genoReadsMapT.begin()->first;
                    if(genoReadsMapT.begin()->second < genoReadsMapT.rbegin()->second) {
                        
                       double ratio = (double)(genoReadsMapT.begin()->second) / genoReadsMapT.rbegin()->second;
                        
                        if (locVarIt->repuitAllLen == 0 || locVarIt->repuitAllLen == locVarIt->repuit.length()) {
                            if (diff <= locVarIt->repuit.length()) {
                                if (ratio < mOptions->mLocVars.locVarOptions.hlRatio1) {
                                    genoReadsMapT.clear();
                                    nogap = true;
                                }
                            } else if (diff <= 2 * locVarIt->repuit.length()) {
                                if (ratio < mOptions->mLocVars.locVarOptions.hlRatio2) {
                                    genoReadsMapT.clear();
                                    nogap = true;
                                }
                            } else {
                                if (ratio < mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                    genoReadsMapT.clear();
                                    nogap = true;
                                }
                            } 
                        } else {
                            if (diff <= locVarIt->repuit.length() || diff <= locVarIt->repuit2.length() || diff <= locVarIt->repuitAllLen) {
                                if (ratio < mOptions->mLocVars.locVarOptions.hlRatio1) {
                                    genoReadsMapT.clear();
                                    nogap = true;
                                }
                            } else if (diff <= 2 * locVarIt->repuit.length() || diff <= 2 * locVarIt->repuit2.length() || diff <= 2 * locVarIt->repuitAllLen) {
                                if (ratio < mOptions->mLocVars.locVarOptions.hlRatio2) {
                                    genoReadsMapT.clear();
                                    nogap = true;
                                }
                            } else {
                                if (ratio < mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                    genoReadsMapT.clear();
                                    nogap = true;
                                }
                            }
                        } 
                                                
                    } else {
                        
                    }

                } else {
                    nogap = true;
                }

                if (nogap) {
                    //get the two trend;
                    int repul = locVarIt->repuit.mStr.length();
                    auto twoPeaksMap = get2Peaks(genoReadsMap, mOptions->mLocVars.locVarOptions.hlRatio1, mOptions->mLocVars.locVarOptions.hlRatio2);

                    if (twoPeaksMap.size() == 1) {
                        genoReadsMapT = twoPeaksMap;
                        twoPeaksMap.clear();
                    } else {

                        double ratio = (double)std::min(twoPeaksMap.begin()->second, twoPeaksMap.rbegin()->second) / std::max(twoPeaksMap.begin()->second, twoPeaksMap.rbegin()->second);

                        int diff = twoPeaksMap.rbegin()->first - twoPeaksMap.begin()->first;
                        if (locVarIt->repuitAllLen == 0 || locVarIt->repuitAllLen == locVarIt->repuit.length()) {
                            if (diff <= locVarIt->repuit.length()) {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            } else if (diff <= 2 * locVarIt->repuit.length()) {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            } else if (diff <= 3 * locVarIt->repuit.length()) {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            } else {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            }
                        } else {
                            if (diff <= locVarIt->repuit.length() || diff <= locVarIt->repuit2.length() || diff <= locVarIt->repuitAllLen) {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            } else if (diff <= 2 * locVarIt->repuit.length() || diff <= 2 * locVarIt->repuit2.length() || diff <= 2 * locVarIt->repuitAllLen) {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            } else if (diff <= 3 * locVarIt->repuit.length() || diff <= 3 * locVarIt->repuit2.length() || diff <= 3 * locVarIt->repuitAllLen) {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            } else {
                                if (ratio >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                                    genoReadsMapT = twoPeaksMap;
                                }
                            }
                        }

                        if (genoReadsMapT.empty()) {
                            auto maxGenoReads = getMaxKeyValue(genoReadsMap); //key max value pair;
                            genoReadsMapT[maxGenoReads.first] = maxGenoReads.second;

                            int occ = 0;
                            for (const auto & it2 : genoReadsMap) {
                                if (it2.second == maxGenoReads.second) {
                                    occ++;
                                }
                            }

                            if (occ == 1) {
                                for (const auto & it2 : genoReadsMap) {
                                    if (it2.first > maxGenoReads.first) {
                                        genoReadsMap2[it2.first] = it2.second;
                                    }
                                }
                                if (!genoReadsMap2.empty()) {
                                    auto maxGenoReads2 = getMaxKeyValue(genoReadsMap2, true);

                                    if (maxGenoReads2.first - maxGenoReads.first <= repul) {
                                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio1) {
                                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                                        }
                                    } else if (maxGenoReads2.first - maxGenoReads.first <= repul * 2) {
                                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio2) {
                                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                                        }
                                    } else if (maxGenoReads2.first - maxGenoReads.first <= repul * 3) {
                                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio2 / 2) {
                                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                                        }
                                    } else {
                                        if ((double) maxGenoReads2.second / maxGenoReads.second >= mOptions->mLocVars.locVarOptions.hlRatio2 / 4) {
                                            genoReadsMapT[maxGenoReads2.first] = maxGenoReads2.second;
                                        }
                                    }
                                }

                                genoReadsMap2.clear();

                            } else if (occ > 1) {
                                for (const auto & it2 : genoReadsMap) {
                                    if (it2.second == maxGenoReads.second && it2.first != maxGenoReads.first) {
                                        genoReadsMapT[it2.first] = maxGenoReads.second;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        genoReadsMap.clear();

        for (auto & it2 : it.second) {
            auto tmpVar = genoReadsMapT.find(it2.second.baseLocVar.effectiveLen);
            it2.second.baseLocVar.totalReads = (tmpVar == genoReadsMapT.end() ? 0 : tmpVar->second);

            *fout << basename(mOptions->prefix) << "\t" << it.first << "\t" << it2.second.baseLocVar.repuitAll.mStr << "\t" <<
                    it2.second.baseLocVar.mraBase << "\t" <<
                    it2.second.baseLocVar.mraName << "\t" << it2.second.baseLocVar.mra.mStr.length() << "\t" <<
                    it2.second.baseLocVar.effectiveLen << "\t" << it2.second.numReads << "\t" <<
                    (it2.second.baseLocVar.totalReads > 0 ? "Y" : "N") << "\t" <<
                    it2.second.baseLocVar.ff.mStr << "\t" << it2.second.baseLocVar.mra.mStr << "\t"
                    << it2.second.baseLocVar.rf.mStr << "\t";
            if (!locVarIt->refSnpsSetffMap[basename(mOptions->prefix)].empty()) {
                int count = 0;
                for (const auto & snps : locVarIt->refSnpsSetffMap[basename(mOptions->prefix)]) {
                    if (locVarIt->ff.mStr[snps] != it2.second.baseLocVar.ff.mStr[snps]) {
                        *fout << snps << locVarIt->ff.mStr[snps] << "|" << it2.second.baseLocVar.ff.mStr[snps] << ";";
                        count++;
                    }
                }
                if (count == 0) {
                    *fout << "NA;";
                }
            } else {
                *fout << "NA;";
            }
            *fout << "\t";

            if (!locVarIt->refSnpsSetrfMap[basename(mOptions->prefix)].empty()) {
                int count = 0;
                for (const auto & snps : locVarIt->refSnpsSetrfMap[basename(mOptions->prefix)]) {
                    if (locVarIt->rf.mStr[snps] != it2.second.baseLocVar.rf.mStr[snps]) {
                        *fout << snps << locVarIt->rf.mStr[snps] << "|" << it2.second.baseLocVar.rf.mStr[snps] << ";";
                        count++;
                    }
                }
                if (count == 0) {
                    *fout << "NA;";
                }
            } else {
                *fout << "NA;";
            }
            *fout << "\n";
        }
        genoReadsMapT.clear();
    }

    fout->flush();
    fout->clear();
    fout->close();
    if (fout) {
        delete fout;
        fout = NULL;
    }
    if (mOptions->verbose) loginfo("End writing genotype table!");
    
    return sortedAllGenotypeMapVec;
}

std::vector<std::pair<std::string, Genotype>> SsrScanner::sortGenotypeMap(std::map<std::string, Genotype> & genoMap) {
    std::vector<std::pair < std::string, Genotype>> sVec;

    std::copy(genoMap.begin(), genoMap.end(),
            std::back_inserter<std::vector<std::pair < std::string, Genotype>>>(sVec));

    std::sort(sVec.begin(), sVec.end(),
            [](const std::pair<std::string, Genotype> & l, const std::pair<std::string, Genotype> & r) {
                if (l.second.baseLocVar.effectiveLen != r.second.baseLocVar.effectiveLen) {
                    return l.second.baseLocVar.effectiveLen < r.second.baseLocVar.effectiveLen;
                }
                return l.first.length() < r.first.length();
            });

    return sVec;
}

std::size_t SsrScanner::mutationMatchFR(std::string target, std::string query, int mis, bool rev) {//target is read, query is flanking region
    std::size_t pos = 0;
    int misMtch = 0;
    if (rev) {

        if (target.length() > query.length()) {

            if (edit_distance(target.substr(target.length() - query.length()), query) <= (uint32) mis) {
                return static_cast<std::size_t> (target.length() - query.length());
            } else {
                int p = 0;
                for (int i = query.length() - 1; i > 0; i--) {
                    if (query[i] != target[target.length() - 1 - p]) {
                        query[i] = target[target.length() - 1 - p];
                        misMtch++;
                        pos = target.rfind(query);
                        if (pos == std::string::npos) {
                            if (misMtch > mis || misMtch > query.length()) {
                                return 0;
                            }
                        } else {
                            if (misMtch <= mis) {
                                return pos;
                            } else {
                                return 0;
                            }
                        }
                    }
                    p++;
                }
                return 0;
            }

        } else {
            int p = 0;
            for (int i = query.length() - 1; i > 0; i--) {
                if (query[i] != target[target.length() - 1 - p]) {
                    query[i] = target[target.length() - 1 - p];
                    misMtch++;
                    pos = target.rfind(query);
                    if (pos == std::string::npos) {
                        if (misMtch > mis || misMtch > query.length()) {
                            return 0;
                        }
                    } else {
                        if (misMtch <= mis) {
                            return pos;
                        } else {
                            return 0;
                        }
                    }
                }
                p++;
            }
            return 0;
        }
    } else {
        if (target.length() > query.length()) {

            if (edit_distance(target.substr(0, query.length()), query) <= (uint32) mis) {
                return 0;
            } else {
                for (int i = 0; i != query.length(); i++) {
                    if (query[i] != target[i]) {
                        query[i] = target[i];
                        misMtch++;
                        pos = target.find(query);
                        if (pos == std::string::npos) {
                            if (misMtch > mis || misMtch > query.length()) {
                                return 1;
                            }
                        } else {
                            if (misMtch <= mis) {
                                return 0;
                            } else {
                                return 1;
                            }
                        }
                    }
                }
                return 1;
            }

        } else {
            for (int i = 0; i != query.length(); i++) {
                if (query[i] != target[i]) {
                    query[i] = target[i];
                    misMtch++;
                    pos = target.find(query);
                    if (pos == std::string::npos) {
                        if (misMtch > mis || misMtch > query.length()) {
                            return 1;
                        }
                    } else {
                        if (misMtch <= mis) {
                            return 0;
                        } else {
                            return 1;
                        }
                    }
                }
            }
            return 1;
        }
    }
}

void SsrScanner::preAnalyze(Read* & r1, std::size_t & ffpos, std::size_t & rfpos, bool & mraAnalyze) {
    int misF = locVarIt->ff.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR == 0 ? 1 : (locVarIt->ff.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR);
    int misR = locVarIt->rf.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR == 0 ? 1 : (locVarIt->rf.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR);

    if (ffpos != std::string::npos && rfpos != std::string::npos) {
        if (static_cast<int> (ffpos) + locVarIt->ff.mStr.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) <= static_cast<int> (rfpos)) {
            std::pair<size_t, int> tmp = analyzeMRA(r1->mSeq.mStr, locVarIt->repuit.mStr, ffpos, rfpos);
            if (tmp.second >= mOptions->mLocVars.locVarOptions.minNSSRUnit * locVarIt->repuit.mStr.length()) {
                ffpos += static_cast<std::size_t> (locVarIt->ff.mStr.length());
            } else {
                mraAnalyze = true;
            }
        } else {
            ffpos = 0;
            rfpos = static_cast<std::size_t> (r1->mSeq.length() - 1);
            mraAnalyze = true;
        }
    } else if (ffpos == std::string::npos && rfpos != std::string::npos) {
        if (locVarIt->ff.mStr.empty()) {
            ffpos = 0;
            if (mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()) > static_cast<int> (rfpos)) {
                rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                mraAnalyze = true;
            }
        } else {
            ffpos = mutationMatchFR(r1->mSeq.mStr, locVarIt->ff.mStr, misF, false);

            if (ffpos == 0) {
                ffpos = static_cast<std::size_t> (locVarIt->ff.mStr.length());
                if (rfpos < ffpos + static_cast<std::size_t> (mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()))) {
                    rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                    mraAnalyze = true;
                }
            } else {
                ffpos = 0;
                if (rfpos < ffpos + static_cast<std::size_t> (mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()))) {
                    rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                }
                mraAnalyze = true;
            }

        }
    } else if (ffpos != std::string::npos && rfpos == std::string::npos) {
        if (locVarIt->rf.mStr.empty()) {
            rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
        } else {

            ffpos += static_cast<std::size_t> (locVarIt->ff.mStr.length());

            rfpos = mutationMatchFR(r1->mSeq.mStr, locVarIt->rf.mStr, misR, true);

            if (rfpos == 0) {
                rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                mraAnalyze = true;
            } else {
                if (rfpos < ffpos + static_cast<std::size_t> (mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()))) {
                    rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                    mraAnalyze = true;
                }
            }
        }

    } else if (ffpos == std::string::npos && rfpos == std::string::npos) {

        if (locVarIt->ff.mStr.empty() && locVarIt->rf.mStr.empty()) {
            ffpos = 0;
            rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
        } else if (!locVarIt->ff.mStr.empty() && locVarIt->rf.mStr.empty()) {
            rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);

            ffpos = mutationMatchFR(r1->mSeq.mStr, locVarIt->ff.mStr, misF, false);

            if (ffpos == 0) {
                ffpos = static_cast<std::size_t> (locVarIt->ff.mStr.length());
            } else {
                ffpos = 0;
                mraAnalyze = true;
            }

        } else if (locVarIt->ff.mStr.empty() && !locVarIt->rf.mStr.empty()) {
            ffpos = 0;

            rfpos = mutationMatchFR(r1->mSeq.mStr, locVarIt->rf.mStr, misR, true);

            if (rfpos == 0) {
                rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                mraAnalyze = true;
            } else {
                if (rfpos < ffpos + static_cast<std::size_t> (mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()))) {
                    rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                    mraAnalyze = true;
                }
            }

        } else if (!locVarIt->ff.mStr.empty() && !locVarIt->rf.mStr.empty()) {
            ffpos = mutationMatchFR(r1->mSeq.mStr, locVarIt->ff.mStr, misF, false);
            rfpos = mutationMatchFR(r1->mSeq.mStr, locVarIt->rf.mStr, misR, true);

            if (ffpos == 0 && rfpos == 0) {
                ffpos = static_cast<std::size_t> (locVarIt->ff.mStr.length());
                rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                mraAnalyze = true;
            } else if (ffpos != 0 && rfpos == 0) {
                ffpos = 0;
                rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                mraAnalyze = true;
            } else if (ffpos == 0 && rfpos != 0) {
                ffpos = static_cast<std::size_t> (locVarIt->ff.mStr.length());
                if (rfpos < static_cast<std::size_t> (locVarIt->ff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length()))) {
                    rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                    mraAnalyze = true;
                }
            } else {
                ffpos = 0;
                if (static_cast<std::size_t> (locVarIt->ff.length() + mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length())) > rfpos) {
                    rfpos = static_cast<std::size_t> (r1->mSeq.mStr.length() - 1);
                }
                mraAnalyze = true;
            }
        }
    }
}

std::string SsrScanner::scanVar (Read* & r1) {
    ss.str("");
    readSeq = r1->mSeq.mStr.c_str();
    readLength = r1->mSeq.length();
    readName = r1->mName;
    returnedlocus.clear();
    checkLoci = true;
    
    std::map<std::string, std::tuple<int, int, int>> locMap; //loci, seq score, trimposF, trimposR;

    for (auto & it : mOptions->mLocVars.refLocMap) {
        if (r1->mSeq.length() > (it.second.fp.length() + it.second.rp.length())) {
            bool goRP = false;
            int trimPosF = 0;
            int fpMismatches = (int) edit_distance(it.second.fp.mStr, r1->mSeq.mStr.substr(0, it.second.fp.length()));

            fpData = it.second.fp.mStr.c_str();
            fpLength = it.second.fp.length();
            auto endBoolF = doPrimerAlignment(fpData, fpLength, it.second.name, readSeq, readLength, r1->mName, true);

            if(std::min(fpMismatches, get<0>(endBoolF)) > mOptions->mLocVars.locVarOptions.maxMismatchesPSeq) continue;

            if(get<0>(endBoolF) <= fpMismatches){
              fpMismatches = get<0>(endBoolF);
              if ((get<1>(endBoolF) + it.second.ff.length() + it.second.rf.length() + it.second.rp.length()) <= r1->length()) {
                trimPosF = get<1>(endBoolF);
                goRP = true;
               }
            } else {
                trimPosF = it.second.fp.length();
                goRP = true;
            }

            if (goRP && (r1->length() >= it.second.rp.length() + it.second.rf.length() + it.second.ff.length())) {
                int rpMismatches = (int) edit_distance(it.second.rp.mStr, r1->mSeq.mStr.substr(r1->mSeq.length() - it.second.rp.length()));
                rpData = it.second.rp.mStr.c_str();
                rpLength = it.second.rp.length();
                auto endBoolR = doPrimerAlignment(rpData, rpLength, it.second.name, readSeq, readLength, r1->mName, true);

                if(std::min(rpMismatches, get<0>(endBoolR)) > mOptions->mLocVars.locVarOptions.maxMismatchesPSeq) continue;

                if(get<0>(endBoolR) <= rpMismatches){
                    rpMismatches = get<0>(endBoolR);
                    if (get<2>(endBoolR) && (get<1>(endBoolR) <= readLength)) {
                        //if(bprint) nnumberR++;
                        locMap[it.second.name] = std::make_tuple((fpMismatches + get<0>(endBoolR)), trimPosF, (get<1>(endBoolR) - trimPosF - it.second.rp.mStr.length()));
                        //locMap[it.second.name] = std::make_pair((fpMismatches + get<0>(endBool)), get<1>(endBool) - it.second.fp.mStr.length() - it.second.rp.mStr.length());
                    }
                } else {
                     locMap[it.second.name] = std::make_tuple((fpMismatches + rpMismatches), trimPosF, (readLength - trimPosF - it.second.rp.mStr.length()));
                    //locMap[it.second.name] = std::make_pair((fpMismatches + rpMismatches), readLength - it.second.fp.mStr.length() - it.second.rp.mStr.length());
                    //should punish if there are mismatches in fp; 
                    //and could not be the best if there are indel in the head of the rp
                } 
            }
        }
    }

    if (locMap.empty()) {
        returnedlocus = "_failed";
    } else {
        std::string locName = "";
        if (locMap.size() == 1) {
            locName = locMap.begin()->first;
        } else {
            std::vector<int> seqScoreVec;
            for (const auto & it : locMap) {
                seqScoreVec.push_back(get<0>(it.second));
            }
            auto minValue = *std::min_element(seqScoreVec.begin(), seqScoreVec.end());
            //warning, what if there are multiple identical values
            seqScoreVec.clear();

            for (const auto & it : locMap) {
                if (get<0>(it.second) == minValue) {
                    locName = it.first;
                    ss.str();
                    break; //warning, what if there are multiple identical values
                }
            }
        }

        locVarIt = &(mOptions->mLocVars.refLocMap[locName]);
        r1->trimFront(get<1>(locMap[locName]));
        r1->resize(get<2>(locMap[locName]));
        locMap.clear();

        minReadLen = std::max(locVarIt->ff.length(), locVarIt->rf.length()) +
                mOptions->mLocVars.locVarOptions.minNSSRUnit * (locVarIt->repuit.mStr.length() + locVarIt->repuit2.mStr.length());
        if (r1->length() < minReadLen) {
            returnedlocus = locVarIt->name + "_failed";
            return (returnedlocus);
        }

        std::map<std::string, std::map < std::string, Genotype>>::iterator itGenotypeMap = tmpAllGenotypeMap.find(locVarIt->name);
        std::map<std::string, Genotype> tmpGenotypeMap;
        Genotype* tmpGenotype = new Genotype();
        if (itGenotypeMap == tmpAllGenotypeMap.end()) {
            if (r1->mSeq.mStr == locVarIt->effectiveSeq.mStr) {//identical seq with ref
                tmpGenotype->baseLocVar = *locVarIt;
                tmpGenotype->numReads++;
                tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
                tmpAllGenotypeMap[locVarIt->name] = tmpGenotypeMap;
                returnedlocus = locVarIt->name + "_true";
                if (mOptions->mEdOptions.printRes) {
                    //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                }
            } else {
                bool mraAnalyze = false;
                std::size_t ffpos = r1->mSeq.mStr.find(locVarIt->ff.mStr);
                std::size_t rfpos = r1->mSeq.mStr.rfind(locVarIt->rf.mStr);
                preAnalyze(r1, ffpos, rfpos, mraAnalyze);
                if (mraAnalyze) {
                    //return the pos and the length of mra after the pos
                    std::pair<size_t, int> mraPosLenReadF = analyzeMRA(r1->mSeq.mStr, locVarIt->repuit.mStr, ffpos, rfpos);

                    ss << "analyzeMRA: " << mraPosLenReadF.first << " : " << mraPosLenReadF.second << " -> " << r1->mSeq.mStr.substr(mraPosLenReadF.first, mraPosLenReadF.second) << "\n";

                    if (mraPosLenReadF.second > 0) {
                        tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, mraPosLenReadF, locVarIt->repuit.mStr);
                        tmpGenotype->numReads++;
                        tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
                        tmpAllGenotypeMap[locVarIt->name] = tmpGenotypeMap;
                        if (mOptions->mEdOptions.printRes) {
                            //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                        }
                        returnedlocus = locVarIt->name + "_true";
                    } else {
                        returnedlocus = locVarIt->name + "_failed";
                    }
                } else {
                    tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, ffpos, rfpos, locVarIt->repuit.mStr);
                    tmpGenotype->numReads++;
                    tmpGenotypeMap[r1->mSeq.mStr] = *tmpGenotype;
                    tmpAllGenotypeMap[locVarIt->name] = tmpGenotypeMap;
                    returnedlocus = locVarIt->name + "_true";
                }

                if (mOptions->mEdOptions.printRes) {
                    // cCout(ss.str(), 'y');
                }
            }

        } else {
            std::map<std::string, Genotype>::iterator tmpGenotypeMapIT = itGenotypeMap->second.find(r1->mSeq.mStr);

            if (tmpGenotypeMapIT == itGenotypeMap->second.end()) {
                if (r1->mSeq.mStr == locVarIt->effectiveSeq.mStr) {//identical seq with ref
                    tmpGenotype->baseLocVar = *locVarIt;
                    tmpGenotype->numReads++;
                    itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
                    if (mOptions->mEdOptions.printRes) {
                        //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                    }
                    returnedlocus = locVarIt->name + "_true";
                } else {
                    bool mraAnalyze = false;
                    std::size_t ffpos = r1->mSeq.mStr.find(locVarIt->ff.mStr);
                    std::size_t rfpos = r1->mSeq.mStr.rfind(locVarIt->rf.mStr);
                    preAnalyze(r1, ffpos, rfpos, mraAnalyze);
                    if (mraAnalyze) {
                        //return the pos and the length of mra after the pos
                        std::pair<size_t, int> mraPosLenReadF = analyzeMRA(r1->mSeq.mStr, locVarIt->repuit.mStr, ffpos, rfpos);
                        if (mraPosLenReadF.second > 0) {
                            tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, mraPosLenReadF, locVarIt->repuit.mStr);
                            tmpGenotype->numReads++;
                            itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
                            if (mOptions->mEdOptions.printRes) {
                                //ss << "Genotype: " << tmpGenotype.numReads << " : " << tmpGenotype.baseLocVar.mra.mStr << "\n";
                            }
                            returnedlocus = locVarIt->name + "_true";
                        } else {
                            returnedlocus = locVarIt->name + "_failed";
                        }
                    } else {
                        tmpGenotype->baseLocVar = LocVar(r1->mSeq.mStr, ffpos, rfpos, locVarIt->repuit.mStr);
                        tmpGenotype->numReads++;
                        itGenotypeMap->second[r1->mSeq.mStr] = *tmpGenotype;
                        returnedlocus = locVarIt->name + "_true";
                    }
                    if (mOptions->mEdOptions.printRes) {
                        //cCout(ss.str(), 'y');
                    }
                }
            } else {
                tmpGenotypeMapIT->second.numReads++;
                returnedlocus = locVarIt->name + "_true";
            }
        }


        if (tmpGenotype) {
            delete tmpGenotype;
            tmpGenotype = nullptr;
        }
    }
    //}
    ss.str();
    return returnedlocus;
}

bool SsrScanner::simpleAnalyzeFlankingRegion(Read* & r1) {

    if (locVarIt->ff.length() == 0 && locVarIt->rf.length() == 0) {
        return true;
    } else if (locVarIt->ff.length() != 0 && locVarIt->rf.length() == 0) {
        int misF = locVarIt->ff.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR == 0 ? 1 : (locVarIt->ff.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR);
        unsigned int edf = edit_distance(locVarIt->ff.mStr, r1->mSeq.mStr.substr(0, locVarIt->ff.length()));
        if (misF >= ((int) edf)) {
            return true;
        } else {
            return false;
        }
    } else if (locVarIt->ff.length() == 0 && locVarIt->rf.length() != 0) {
        int misR = locVarIt->rf.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR == 0 ? 1 : (locVarIt->rf.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR);
        unsigned int edr = edit_distance(locVarIt->rf.mStr, r1->mSeq.mStr.substr(r1->mSeq.length() - locVarIt->rf.length()));
        if (misR >= ((int) edr)) {
            return true;
        } else {
            return false;
        }
    } else {
        int misF = locVarIt->ff.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR == 0 ? 1 : (locVarIt->ff.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR);
        int misR = locVarIt->rf.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR == 0 ? 1 : (locVarIt->rf.mStr.length() * mOptions->mLocVars.locVarOptions.maxMismatchesPer4FR);
        unsigned int edf = edit_distance(locVarIt->ff.mStr, r1->mSeq.mStr.substr(0, locVarIt->ff.length()));
        unsigned int edr = edit_distance(locVarIt->rf.mStr, r1->mSeq.mStr.substr((r1->mSeq.length() - locVarIt->rf.length())));

        if (misF >= ((int) edf) && misR >= ((int) edr)) {
            return true;
        } else {
            return false;
        }
    }

}

std::tuple<int, int, bool> SsrScanner::doPrimerAlignment(const char* & qData, int qLength, const std::string & qName,
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

void SsrScanner::printAlignment2(const char* query, const char* target,
        const unsigned char* alignment, const int alignmentLength,
        const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx = -1;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                printf("-");
            else
                printf("%c", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);

        // match / mismatch
        printf("   ");
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");

        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                printf("-");
            else
                printf("%c", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
    }
}
