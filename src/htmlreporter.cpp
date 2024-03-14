#include "htmlreporter.h"
#include <chrono>
#include <memory.h>
#include <valarray>

extern string command;

HtmlReporter::HtmlReporter(Options* opt) {
    mOptions = opt;
    mDupHist = NULL;
    mDupRate = 0.0;
}

HtmlReporter::~HtmlReporter() {
}

void HtmlReporter::setDupHist(int* dupHist, double* dupMeanGC, double dupRate) {
    mDupHist = dupHist;
    mDupMeanGC = dupMeanGC;
    mDupRate = dupRate;
}

void HtmlReporter::setInsertHist(long* insertHist, int insertSizePeak) {
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

void HtmlReporter::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, std::string & marker, std::vector<std::pair<std::string, Genotype>> &outGenotype, Options*& mOptions) {
    for (auto & it : outGenotype) {
        ofs << "<tr>";
        ofs << "<td>" + marker + "</td>" +
                "<td>" + it.second.baseLocVar.repuitAll.mStr + "</td>" +
                "<td>" + it.second.baseLocVar.mraBase + "</td>" +
                "<td>" + std::to_string(it.second.baseLocVar.mra.mStr.length()) + "</td>" +
                "<td bgcolor=" + (it.second.baseLocVar.totalReads > 0 ? "'green'>" : "'transparent'>") +
                std::to_string(it.second.baseLocVar.effectiveLen) + "</td>" +
                "<td><font color='" + (it.second.numReads < mOptions->mLocVars.locVarOptions.minWarningSeqs ? "red" : "black") + "'>" +
                std::to_string(it.second.numReads) + "</font></td>" +
                //"<td>" + std::to_string(it.second.numReads) + "</td>" +
                "<td align='right'>" + highligher(it.second.baseLocVar.ff.mStr, it.second.baseLocVar.snpsMapff) + "</td>" +
                "<td align='center'>" + it.second.baseLocVar.mraName + "</td>" +
                "<td align='left'>" + highligher(it.second.baseLocVar.rf.mStr, it.second.baseLocVar.snpsMaprf) + "</td>";
        ofs << "</tr>\n";
    }
}

void HtmlReporter::outputRow(ofstream& ofs, LocSnp2& locSnp, bool align = true, int num = 5) {
    int i = 1;
    if (align) {
        for (int ii = 0; ii < locSnp.seqVarVec.size(); ii++) {
            if(ii > num) break;
            std::string hap = ""; //haplotype str or indel;
            std::string zygosity = ""; // homo, heter, inconclusive, or seq error;
            std::string snpsStr = ""; //snp str;
            std::string haploRatio = "NA";
            std::string fc = "black";
            std::string bc = "transparent";

            if (ii == 0) {
                zygosity = locSnp.genoStr3;
                fc = "white";
                haploRatio = std::to_string(locSnp.ratioHaplo);
                if (locSnp.status.first.first) {
                    hap = "indel";
                    bc = "gray";
                } else {
                    hap = locSnp.getHaploStr();
                    snpsStr = locSnp.getSnpStr();
                    bc = "olive";
                }
            } else if (ii == 1) {
                if (locSnp.genoStr3 == "homo") {
                    zygosity = "seq error";
                    if (!locSnp.seqVarVec.at(ii).indel) {
                        hap = locSnp.getHaploStr(ii);
                        snpsStr = locSnp.getSnpStr(ii);
                    } else {
                        hap = "indel";
                    }
                } else {
                    zygosity = locSnp.genoStr3;
                    fc = "white";
                    haploRatio = std::to_string(1 - locSnp.ratioHaplo);
                    if (locSnp.status.first.second) {
                        hap = "indel";
                        bc = "gray";
                    } else {
                        hap = locSnp.getHaploStr(true);
                        snpsStr = locSnp.getSnpStr(true);
                        bc = "olive";
                    }
                }
            } else {
                zygosity = "seq error";
                if (!locSnp.seqVarVec.at(ii).indel) {
                    hap = locSnp.getHaploStr(ii);
                    snpsStr = locSnp.getSnpStr(ii);
                } else {
                    hap = "indel";
                }
            }

            ofs << "<tr>";
            ofs << "<td>" + std::to_string(i) + "</td>" +
                    "<td>" + locSnp.name + "</td>" +
                    "<td>" + snpsStr + "</td>" +
                    "<td bgcolor='" + bc + "'>" + //haplotype
                    "<font color='" + fc + "'>" + hap + "</font></td>" +
                    "<td bgcolor='" + bc + "'>" + //haplotype
                    "<font color='" + fc + "'>" + std::to_string(locSnp.seqVarVec.at(ii).numReads) + "</font></td>" +
                    "<td bgcolor='" + bc + "'>" + //
                    "<font color='" + fc + "'>" + haploRatio + "</font></td>" + //homo or heter or inconclusive;
                    "<td bgcolor='" + bc + "'>" + //
                    "<font color='" + fc + "'>" + zygosity + "</font></td>" + //homo or heter or inconclusive;
                    "<td>" + (ii == 0 ? locSnp.getRatioStr() : "") + "</td>" +
                    "<td>" + std::to_string(locSnp.getReadsVarPer(ii)) + "</td>" +
                    "<td>" + std::to_string(locSnp.totReads) + "</td>" +
                    "<td>" + std::to_string(locSnp.seqVarVec.at(ii).seq.length()) + "</td>" +
                    "<td align='center'>" + highligher(locSnp, false, locSnp.ref.mStr, locSnp.seqVarVec.at(ii).seq, locSnp.seqVarVec.at(ii).snpSet) + 
                    "</td>";
            ofs << "</tr>\n";
            i++;
        }
    } else {
        for (const auto & it : locSnp.ssnpsMap) {
            //std::string fc = it.second.color == 'o' ? "black" : "white";
            std::string fc = "white";
            std::string bc = it.second.color == 'g' ? "green" : (it.second.color == 'r' ? "red" : "orange");
            ofs << "<tr>";
            ofs << "<td>" + std::to_string(i) + "</td>" + //ID
                    "<td>" + std::to_string(it.first + locSnp.ft.length()) + "</td>" + //Position
                    "<td bgcolor='" + bc + "'>" + //Genotype
                    "<font color='" + fc + "'>" + it.second.snp1 + "|" + it.second.snp2 + "</font></td>" +
                    "<td bgcolor='" + bc + "'>" + //Putative Genotype
                    "<font color='" + fc + "'>" + (locSnp.genoStr3 != "inconclusive" ? "yes" : "inconclusive") + "</font></td>" +
                    "<td>" + std::to_string(locSnp.getHaploReads()) + "</td>" +
                    "<td>" + std::to_string(locSnp.getHaploReads(true)) + "</td>" +
                    "<td>" + std::to_string(locSnp.ratioHaplo) + "</td>" +
                    "<td>" + std::to_string(locSnp.totHaploReads) + "</font></td>";
            ofs << "</tr>\n";
            i++;
        }
    }
}

//void HtmlReporter::outputRow(ofstream& ofs, LocSnp2& locSnp, bool align = true) {
//    int i = 1;
//    if (align) {
//        for (auto & it : locSnp.genoMap) {
//            std::string hap = "";
//            std::string fc = "black";
//            std::string bc = "transparent";
//       
//            if(it.second.genoStr8 == "seqerr"){
//                 hap =  "seq error";
//            } else {
//                if (it.second.genoStr8 == "indel1") {
//                    hap = "indel1";
//                    bc = "gray";
//                } else if (it.second.genoStr8 == "indel2") {
//                    hap = "indel2";
//                    bc = "gray";
//                } else {
//                    hap = it.second.haploStr;
//                    if (it.second.genoStr8 == "inHeter1") {
//                        bc = "gray";
//                    } else if (it.second.genoStr8 == "inHeter2") {
//                        bc = "gray";
//                    } else {
//                        bc = "olive";
//                    }
//                }
//                fc = "white";
//            }
//            ofs << "<tr>";
//            ofs << "<td>" + std::to_string(i) + "</td>" + 
//                    "<td>" + locSnp.name + "</td>" +
//                    "<td>" + it.second.snpsStr + "</td>" +
//                    "<td bgcolor='" + bc + "'>" + //Genotype
//                    "<font color='" + fc + "'>" + hap + "</font></td>" +
//                    "<td>" + std::to_string(it.second.numReads) + "</td>" +
//                    "<td>" + std::to_string((double) it.second.numReads * 100 / locSnp.totReads) + "</td>" +
//                    "<td>" + std::to_string(locSnp.totReads) + "</td>" +
//                    "<td>" + std::to_string(it.first.length()) + "</td>" +
//                    "<td align='center'>" + highligher(locSnp, false, locSnp.ref.mStr, it.first, it.second.snpPosSet) + "</td>";
//            ofs << "</tr>\n";
//            i++;
//        }
//    } else {
//        for (const auto & it : locSnp.snpsMap) {
//            std::string fc = it.second.color == "transparent" ? "black" : "white";
//            ofs << "<tr>";
//            ofs << "<td>" + std::to_string(i) + "</td>" + //ID
//                    "<td>" + std::to_string(it.first + locSnp.ft.length()) + "</td>" + //Position
//                    "<td bgcolor='" + (it.second.color) + "'>" + //Genotype
//                    "<font color='" + fc + "'>" + it.second.snp1 + "|" + it.second.snp2 + "</font></td>" +
//                    "<td bgcolor='" + (it.second.color) + "'>" + //Putative Genotype
//                    "<font color='" + fc + "'>" + (it.second.genoStr3 != "inconclusive" ? "yes" : "inconclusive") + "</font></td>" +
//                    "<td>" + std::to_string(it.second.reads1) + "</td>" +
//                    "<td>" + std::to_string(it.second.reads2) + "</td>" +
//                    "<td>" + std::to_string(it.second.ratio) + "</td>" +
//                    "<td>" + std::to_string(locSnp.totHaploReads) + "</font></td>";
//            ofs << "</tr>\n";
//            i++;
//        }
//    }
//}

string HtmlReporter::formatNumber(long number) {
    double num = (double) number;
    string unit[6] = {"", "K", "M", "G", "T", "P"};
    int order = 0;
    while (num > 1000.0) {
        order += 1;
        num /= 1000.0;
    }

    if (order == 0)
        return to_string(number);
    else
        return to_string(num) + " " + unit[order];
}

string HtmlReporter::getPercents(long numerator, long denominator) {
    if (denominator == 0)
        return "0.0";
    else
        return to_string((double) numerator * 100.0 / (double) denominator);
}

void HtmlReporter::printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    long pre_total_reads = preStats1->getReads();
    if (preStats2)
        pre_total_reads += preStats2->getReads();

    long pre_total_bases = preStats1->getBases();
    if (preStats2)
        pre_total_bases += preStats2->getBases();

    long pre_q20_bases = preStats1->getQ20();
    if (preStats2)
        pre_q20_bases += preStats2->getQ20();

    long pre_q30_bases = preStats1->getQ30();
    if (preStats2)
        pre_q30_bases += preStats2->getQ30();

    long pre_total_gc = preStats1->getGCNumber();
    if (preStats2)
        pre_total_gc += preStats2->getGCNumber();

    long post_total_reads = postStats1->getReads();
    if (postStats2)
        post_total_reads += postStats2->getReads();

    long post_total_bases = postStats1->getBases();
    if (postStats2)
        post_total_bases += postStats2->getBases();

    long post_q20_bases = postStats1->getQ20();
    if (postStats2)
        post_q20_bases += postStats2->getQ20();

    long post_q30_bases = postStats1->getQ30();
    if (postStats2)
        post_q30_bases += postStats2->getQ30();

    long post_total_gc = postStats1->getGCNumber();
    if (postStats2)
        post_total_gc += postStats2->getGCNumber();

    string sequencingInfo = mOptions->isPaired() ? "paired end" : "single end";
    if (mOptions->isPaired()) {
        sequencingInfo += " (" + to_string(preStats1->getCycles()) + " cycles + " + to_string(preStats2->getCycles()) + " cycles)";
    } else {
        sequencingInfo += " (" + to_string(preStats1->getCycles()) + " cycles)";
    }

    ofs << endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Data QC summary <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='summary' style='display:none'>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n";
    ofs << "<div id='general'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "sequencing:", sequencingInfo);

    // report read length change
    if (mOptions->isPaired()) {
        outputRow(ofs, "mean length before filtering:", to_string(preStats1->getMeanLength()) + "bp, " + to_string(preStats2->getMeanLength()) + "bp");
    } else {
        outputRow(ofs, "mean length before filtering:", to_string(preStats1->getMeanLength()) + "bp");
        outputRow(ofs, "mean length after filtering:", to_string(postStats1->getMeanLength()) + "bp");
    }

    if (mOptions->duplicate.enabled) {
        string dupStr = to_string(mDupRate * 100) + "%";
        if (!mOptions->isPaired())
            dupStr += " (may be overestimated since this is SE data)";
        outputRow(ofs, "duplication rate:", dupStr);
    }
    if (mOptions->isPaired()) {
        outputRow(ofs, "Insert size peak:", mInsertSizePeak);
    }
    if (mOptions->adapterCuttingEnabled()) {
        if (!mOptions->adapter.detectedAdapter1.empty())
            outputRow(ofs, "Detected read1 adapter:", mOptions->adapter.detectedAdapter1);
        if (!mOptions->adapter.detectedAdapter2.empty())
            outputRow(ofs, "Detected read2 adapter:", mOptions->adapter.detectedAdapter2);
    }
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('before_filtering_summary')>Original data</div>\n";
    ofs << "<div id='before_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "total reads:", formatNumber(pre_total_reads));
    outputRow(ofs, "total bases:", formatNumber(pre_total_bases));
    outputRow(ofs, "Q20 bases:", formatNumber(pre_q20_bases) + " (" + getPercents(pre_q20_bases, pre_total_bases) + "%)");
    outputRow(ofs, "Q30 bases:", formatNumber(pre_q30_bases) + " (" + getPercents(pre_q30_bases, pre_total_bases) + "%)");
    outputRow(ofs, "GC content:", getPercents(pre_total_gc, pre_total_bases) + "%");
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('after_filtering_summary')>Clean data used for detection</div>\n";
    ofs << "<div id='after_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "total reads:", formatNumber(post_total_reads));
    outputRow(ofs, "total bases:", formatNumber(post_total_bases));
    outputRow(ofs, "Q20 bases:", formatNumber(post_q20_bases) + " (" + getPercents(post_q20_bases, post_total_bases) + "%)");
    outputRow(ofs, "Q30 bases:", formatNumber(post_q30_bases) + " (" + getPercents(post_q30_bases, post_total_bases) + "%)");
    outputRow(ofs, "GC content:", getPercents(post_total_gc, post_total_bases) + "%");
    ofs << "</table>\n";
    ofs << "</div>\n";

    if (result) {
        ofs << "<div class='subsection_title' onclick=showOrHide('filtering_result')>Filtering result</div>\n";
        ofs << "<div id='filtering_result'>\n";
        result -> reportHtml(ofs, pre_total_reads, pre_total_bases);
        ofs << "</div>\n";
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    if (result && mOptions->adapterCuttingEnabled()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('adapters')><a name='summary'>Adapters <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='adapters' style='display:none'>\n";

        result->reportAdapterHtml(ofs, pre_total_bases);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if (mOptions->duplicate.enabled) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('duplication')><a name='summary'>Duplication <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='duplication' style='display:none'>\n";

        reportDuplication(ofs);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if (mOptions->isPaired()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('insert_size')><a name='summary'>Insert size estimation <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='insert_size' style='display:none'>\n";

        reportInsertSize(ofs, preStats1->getCycles() + preStats2->getCycles() - mOptions->overlapRequire);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

void HtmlReporter::reportInsertSize(ofstream& ofs, int isizeLimit) {
    if (isizeLimit < 1)
        isizeLimit = 1;
    int total = min(mOptions->insertSizeMax, isizeLimit);
    long *x = new long[total];
    double allCount = 0;
    for (int i = 0; i < total; i++) {
        x[i] = i;
        allCount += mInsertHist[i];
    }
    allCount += mInsertHist[mOptions->insertSizeMax];
    double* percents = new double[total];
    memset(percents, 0, sizeof (double)*total);
    if (allCount > 0) {
        for (int i = 0; i < total; i++) {
            percents[i] = (double) mInsertHist[i] * 100.0 / (double) allCount;
        }
    }

    double unknownPercents = (double) mInsertHist[mOptions->insertSizeMax] * 100.0 / (double) allCount;

    ofs << "<div id='insert_size_figure'>\n";
    ofs << "<div class='figure' id='plot_insert_size' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "<div class='sub_section_tips'>This estimation is based on paired-end overlap analysis, and there are ";
    ofs << std::to_string(unknownPercents);
    ofs << "% reads found not overlapped. <br /> The nonoverlapped read pairs may have insert size &lt;" << mOptions->overlapRequire;
    ofs << " or &gt;" << isizeLimit;
    ofs << ", or contain too much sequencing errors to be detected as overlapped.";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'Insert size distribution (" + std::to_string(unknownPercents) + "% reads are with unknown length)', xaxis:{title:'Insert size'}, yaxis:{title:'Read percent (%)'}};\n";
    json_str += "Plotly.newPlot('plot_insert_size', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
}

void HtmlReporter::reportDuplication(ofstream& ofs) {

    ofs << "<div id='duplication_figure'>\n";
    ofs << "<div class='figure' id='plot_duplication' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    int total = mOptions->duplicate.histSize - 2;
    long *x = new long[total];
    double allCount = 0;
    for (int i = 0; i < total; i++) {
        x[i] = i + 1;
        allCount += mDupHist[i + 1];
    }
    double* percents = new double[total];
    memset(percents, 0, sizeof (double)*total);
    if (allCount > 0) {
        for (int i = 0; i < total; i++) {
            percents[i] = (double) mDupHist[i + 1] * 100.0 / (double) allCount;
        }
    }
    int maxGC = total;
    double* gc = new double[total];
    for (int i = 0; i < total; i++) {
        gc[i] = (double) mDupMeanGC[i + 1] * 100.0;
        // GC ratio will be not accurate if no enough reads to average
        if (percents[i] <= 0.05 && maxGC == total)
            maxGC = i;
    }

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Read percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "},";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, maxGC) + "],";
    json_str += "y:[" + Stats::list2string(gc, maxGC) + "],";
    json_str += "name: 'Mean GC ratio (%)  ',";
    json_str += "mode:'lines',";
    json_str += "line:{color:'rgba(255,0,128,1.0)', width:2}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'duplication rate (" + to_string(mDupRate * 100.0) + "%)', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};\n";
    json_str += "Plotly.newPlot('plot_duplication', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
    delete[] gc;
}

void HtmlReporter::reportSex(ofstream & ofs) {

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('sex')><a name='sex'>Sex loci <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='sex'  style='display:none'>\n";

    std::vector<std::string> x_vec{"X", "Y"};
    std::vector<int> y_vec{mOptions->mSex.getHaploVar('x', 0).numReads, mOptions->mSex.getHaploVar('y', 0).numReads};
    std::vector<double> bar_width_vec{0.5, 0.5};

    std::string subsection = "Sex loci: " + mOptions->mSex.sexMarker;
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title = "Sex: " + mOptions->mSex.sexMF + (mOptions->mSex.sexMF == "Male" ? "<br>Y/X = " + std::to_string(mOptions->mSex.YXRatio) : "");

    ofs << "<div class='subsection_title' onclick=showOrHide('" + divName + "')><a name='" + subsection + "'>" + subsection + "<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each allele size will be shown on mouse over.</div>\n";

    ofs << "<div class='figurefull' id='plot_" + divName + "'></div>\n";
    if (!mOptions->mSex.baseErrorMapX.empty() || !mOptions->mSex.baseErrorMapY.empty()) {
        ofs << "<div class='figurehalf' id='plot_sex_e" + divName + "'></div>\n";
    }

    ofs << "<div class='sub_section_tips'>SNPs/artifacts are highlighted in red</div>\n";
    ofs << "<pre overflow: scroll>\n";
    ofs << "<table class='summary_table' style='width:100%'>\n";
    ofs << "<tr style='background:#cccccc'> <td>Sex allele</td><td>N. reads</td><td>Total reads</td><td>Percentage(%)</td><td>Haplotype</td><td>Allele Size</td><td align='center'>Sequence</td></tr>\n";

    if (!mOptions->mSex.sexMarker.empty()) {
        if (!mOptions->mSex.seqVarVecY.empty()) {
            ofs << "<tr style='color:blue'>";
            ofs << "<td>Ref. Y</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>" + std::to_string(mOptions->mSex.refY.length()) + "</td>" <<
                    "<td align='center'>" + highligher(mOptions->mSex.refY.mStr, mOptions->mSex.totSnpSetY) + "</td>";
            ofs << "</tr>";
            
            int ii = 1;
            for (auto it = mOptions->mSex.seqVarVecY.begin(); it != mOptions->mSex.seqVarVecY.end(); ++it) {
                if(ii > mOptions->mLocSnps.mLocSnpOptions.maxRows4Align) break;
                ii++;
                std::string note = "";
                if (it == mOptions->mSex.seqVarVecY.begin()) {
                    if (mOptions->mSex.sexMF == "Male") {
                        note = "haplotype";
                    } else if (mOptions->mSex.sexMF == "Female") {
                        note = "Inconclusive";
                    } else {
                        note = "Inconclusive";
                    }
                } else {
                    note = "seq error";
                }

                ofs << "<tr>";
                ofs << "<td>Y</td>" <<
                        "<td>" + std::to_string(it->numReads) + "</td>" <<
                        "<td>" + std::to_string(mOptions->mSex.totReadsY) + "</td>" <<
                        "<td>" + std::to_string(getPer(it->numReads, mOptions->mSex.totReadsY)) + "</td>" <<
                        "<td>" + note + "</td>" <<
                        "<td>" + std::to_string(it->seq.length()) + "</td>" <<
                        "<td align='center'>" + highligher(it->seq, it->snpSet) + "</td>";
                ofs << "</tr>";
            }
        }

        if (!mOptions->mSex.seqVarVecX.empty()) {
            ofs << "<tr style='color:blue'>";
            ofs << "<td>Ref. X</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>N.A.</td>" <<
                    "<td>" + std::to_string(mOptions->mSex.refX.length()) + "</td>" <<
                    "<td align='center'>" + highligher(mOptions->mSex.refX.mStr, mOptions->mSex.totSnpSetX) + "</td>";
            ofs << "</tr>";
            
            int ii = 1;
            for (auto it = mOptions->mSex.seqVarVecX.begin(); it != mOptions->mSex.seqVarVecX.end(); ++it) {
                if(ii > mOptions->mLocSnps.mLocSnpOptions.maxRows4Align) break;
                ii++;
                std::string note = "";
                if (it == mOptions->mSex.seqVarVecX.begin()) {
                    if (mOptions->mSex.sexMF == "Male") {
                        note = "haplotype";
                    } else if (mOptions->mSex.sexMF == "Female") {
                        note = "haplotype";
                    } else {
                        note = "Inconclusive";
                    }
                } else {
                    if (mOptions->mSex.haplotype) {
                        if (it->seq == mOptions->mSex.getHaploVar('x', 1).seq) {
                            note = "haplotype";
                        } else {
                            note = "seq error";
                        }
                    } else {
                        note = "seq error";
                    }
                }

                ofs << "<tr>";
                ofs << "<td>X</td>" <<
                        "<td>" + std::to_string(it->numReads) + "</td>" <<
                        "<td>" + std::to_string(mOptions->mSex.totReadsX) + "</td>" <<
                        "<td>" + std::to_string(getPer(it->numReads, mOptions->mSex.totReadsX)) + "</td>" <<
                        "<td>" + note + "</td>" <<
                        "<td>" + std::to_string(it->seq.length()) + "</td>" <<
                        "<td align='center'>" + highligher(it->seq, it->snpSet) + "</td>";
                ofs << "</tr>";
            }
        }
    }

    ofs << "</table>\n";
    ofs << "</pre>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;

    string json_str = "var data=[";
    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, x_vec.size()) + "],";
    json_str += "y:[" + Stats::list2string2(y_vec, y_vec.size()) + "],";
    json_str += "type:'bar', textposition: 'auto'";
    json_str += "}";
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "',";
    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(x_vec, x_vec.size()) + "],  title:'" + "Sex Alleles" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    x_vec.clear();
    x_vec.shrink_to_fit();
    y_vec.clear();
    y_vec.shrink_to_fit();
    bar_width_vec.clear();
    bar_width_vec.shrink_to_fit();
    reportSeqError(ofs, divName);

    ofs << "</div>\n";
    ofs << "</div>\n";
}

void HtmlReporter::reportSeqError(ofstream& ofs, std::string & divName) {
    if (mOptions->mSex.baseErrorMapX.empty() && mOptions->mSex.baseErrorMapY.empty()) {
        return;
    }

    ofs << "\n<script type=\"text/javascript\">" << endl;
    std::string json_str = "";

    if (!mOptions->mSex.baseErrorMapY.empty()) {
        std::vector<int> keyV;
        keyV.reserve(mOptions->mSex.baseErrorMapY.size());
        std::vector<double> valueV;
        valueV.reserve(mOptions->mSex.baseErrorMapY.size());
        for (const auto & posIt : mOptions->mSex.baseErrorMapY) {
            keyV.push_back(posIt.first);
            valueV.push_back(posIt.second);
        }

        json_str += "var Ydata = {";
        json_str += "x:[" + Stats::list2string2(keyV, keyV.size()) + "],\n";
        json_str += "y:[" + Stats::list2string2(valueV, valueV.size()) + "],\n";
        json_str += "type: 'scatter', mode: 'markers', marker:{color: 'blue'}, textposition: 'auto', name: 'Y allele'";
        json_str += "};\n";
    }

    if (!mOptions->mSex.baseErrorMapX.empty()) {
        std::vector<int> keyV;
        keyV.reserve(mOptions->mSex.baseErrorMapX.size());
        std::vector<double> valueV;
        valueV.reserve(mOptions->mSex.baseErrorMapX.size());
        for (const auto & posIt : mOptions->mSex.baseErrorMapX) {
            keyV.push_back(posIt.first);
            if (mOptions->mSex.baseErrorMapY.empty()) {
                valueV.push_back(posIt.second);
            } else {
                posIt.second == 0 ? valueV.push_back(posIt.second) : valueV.push_back(-posIt.second);
            }
        }

        json_str += "var Xdata = {";
        json_str += "x:[" + Stats::list2string2(keyV, keyV.size()) + "],\n";
        json_str += "y:[" + Stats::list2string2(valueV, valueV.size()) + "],\n";
        json_str += "type: 'scatter', mode: 'markers', marker:{color: 'red'}, textposition: 'auto', name: 'X allele'";
        json_str += "};\n";
    }
    
    if (mOptions->mSex.baseErrorMapY.empty()) {
        json_str += "var layout = { title: 'Sequence error rate', xaxis:{title:'Position', automargin: true},";
        json_str += "yaxis:{showline:false, zeroline:false, zerolinecolor:'transparent', zerolinewidth:2, title:'Error rate (%)', automargin:true},";
        json_str += "shapes: [";
        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string(mOptions->mSex.aveErrorRateX) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string(mOptions->mSex.aveErrorRateX) + ", ";
        json_str += "line:{color: 'gray', width: 2, dash: 'dash'}";
        json_str += "}]";
        json_str += "};\n";
        json_str += "var data = [Xdata];\n";
    } else {
        json_str += "var layout = { title: 'Sequence error rate', xaxis:{title:'Position', automargin: true},";
        json_str += "yaxis:{showline:false, zeroline:false, zerolinecolor:'transparent', zerolinewidth:2, title:'Error rate (%)', automargin:true},";
        json_str += "shapes: [";
        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string(mOptions->mSex.aveErrorRateY) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string(mOptions->mSex.aveErrorRateY) + ", ";
        json_str += "line:{color: 'gray', width: 2, dash: 'dash'}},\n";
        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string(-mOptions->mSex.aveErrorRateX) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string(-mOptions->mSex.aveErrorRateX) + ", ";
        json_str += "line:{color: 'gray', width: 2, dash: 'dash'}";
        json_str += "}]";
        json_str += "};\n";
        json_str += "var data = [Ydata, Xdata];\n";
    }
    json_str += "Plotly.newPlot('plot_sex_e" + divName + "', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
}

void HtmlReporter::report(std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> &sortedAllGenotypeMapVec,
        FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    ofstream ofs;
    ofs.open(mOptions->htmlFile, ifstream::out);

    printHeader(ofs);

    ofs << "<h1 style='text-align:left;'><a href='https://github.com/seq2sat' target='_blank' style='color:#009900;text-decoration:none;'>Seq2Sat Report</a </h1>" << endl;
    ofs << "<div style='font-size:12px;font-weight:normal;text-align:left;color:#666666;padding:5px;'>" << "Sample: " << basename(mOptions->prefix) << "</div>" << endl;

    if (!mOptions->mSex.sexMarker.empty()) {
        reportSex(ofs);
    }

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('genotype')><a name='genotype'>All genotypes <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='genotype' >\n";

    if (mOptions->mVarType == ssr) {
        reportAllGenotype(ofs, sortedAllGenotypeMapVec);
    } else {
        if(!mOptions->noPlot) reportAllSnps(ofs);
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    printSummary(ofs, result, preStats1, postStats1, preStats2, postStats2);

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('before_filtering')><a name='summary'>Original data <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='before_filtering'  style='display:none'>\n";

    if (preStats1) {
        preStats1 -> reportHtml(ofs, "Original data", "read1");
    }

    if (preStats2) {
        preStats2 -> reportHtml(ofs, "Original data", "read2");
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('after_filtering')><a name='summary'>Clean data used for detection <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='after_filtering'  style='display:none'>\n";

    if (postStats1) {
        string name = "read1";
        postStats1 -> reportHtml(ofs, "Clean data used for detection", name);
    }

    if (postStats2) {
        postStats2 -> reportHtml(ofs, "Clean data used for detection", "read2");
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    printFooter(ofs);

}

void HtmlReporter::reportAllSnps(ofstream& ofs) {
    for (auto & it : mOptions->mLocSnps.refLocMap) {
        if (it.second.seqVarVec.empty()) continue;
        std::string subsection = "Marker: " + it.first;
        std::string divName = replace(subsection, " ", "_");
        divName = replace(divName, ":", "_");
        std::string title = it.first;
        ofs << "<div class='subsection_title' onclick=showOrHide('" + divName + "')><a name='" + subsection + "'>" + subsection + "<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        //ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" + divName + "')>" + subsection + "</a></div>\n";
        ofs << "<div id='" + divName + "'>\n";
        ofs << "<div class='sub_section_tips'>Value of each allele size will be shown on mouse over.</div>\n";
        reportSnpAlignmentTable(ofs, divName, it.second);
        if (!it.second.status.first.second && !it.second.snpPosSetHaplo.empty() && !mOptions->noSnpPlot) {
            reportSnpTablePlot(ofs, divName, it.second);
        }
        ofs << "</div>\n";
    }
}

void HtmlReporter::reportSnpAlignmentTable(ofstream& ofs, std::string & divName, LocSnp2 & locSnp) {
    if (!mOptions->noSnpPlot) {
        ofs << "<div class='figurehalf' id='plot_h" + divName + "'></div>\n";
    }
    
    if (!locSnp.baseErrorMap.empty() && !mOptions->noErrorPlot) {
        ofs << "<div class='figurefull' id='plot_e" + divName + "'></div>\n";
    }

    if (locSnp.seqVarVec.empty()) return;
    ofs << "<div class='sub_section_tips'><font color='red'>Target heter SNPs are highlighted red, </font> <font color='green'>while target homo SNPs are highlighted green, </font> <font color='orange'>new potential SNPs are orange, </font><font color='gray'>sequence artifacts (errors) are gray!</font></div>\n";
    ofs << "<pre overflow: scroll>\n";
    ofs << "<table class='summary_table' style='width:100%'>\n";
    ofs << "<tr style='background:#cccccc'> <td>ID</td><td>Marker</td><td>SNV</td><td>Haplotype</td><td>N. of Reads</td><td>Haplotype Ratio</td><td>Zygosity</td><td>Variant Ratio</td><td>Reads(%)</td><td>Total Reads</td><td>Length</td><td align='center'>Sequence</td></tr>\n";

    ofs << "<tr style='color:blue'>";
    ofs << "<td>0</td>" <<
            "<td>Reference</td>" <<
            "<td>ref</td>" <<
            "<td>N.A.</td>" <<
            "<td>N.A.</td>" <<
            "<td>N.A.</td>" <<
            "<td>N.A.</td>" <<
            "<td>N.A.</td>" <<
            "<td>N.A.</td>" <<
            "<td>N.A.</td>" <<
            "<td>" + std::to_string(locSnp.ft.length() + locSnp.ref.mStr.length() + locSnp.rt.length()) + "</td>" <<
            "<td align='center'>" + highligher(locSnp, true, locSnp.ref.mStr, locSnp.ref.mStr, locSnp.snpPosSet) + "</td>";
    ofs << "</tr>\n";
    outputRow(ofs, locSnp, true, mOptions->mLocSnps.mLocSnpOptions.maxRows4Align);
    ofs << "</table>\n";
    ofs << "</pre>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;

    string json_str = "";
    if (!mOptions->noSnpPlot) {
        // for bar plot of haplotype;
        json_str += "var data=[{";
        json_str += "x:['Allele', 'Allele'],";
        json_str += "y:[" + std::to_string(locSnp.getHaploReads()) + ", " + std::to_string(locSnp.getHaploReads(true)) + "],";
        json_str += "text: ['" + locSnp.getHaploStr() + "', '" + locSnp.getHaploStr(true) + "'],";
        json_str += "width: [0.5, 0.5],";
        json_str += "type:'bar', textposition: 'auto', ";

        if (locSnp.genoStr3 == "homo") {
            json_str += "marker:{color: ['darkgreen', 'darkgreen'], line: {color: 'white', width: 1.5}}";
        } else {
            json_str += "marker:{color: ['goldenrod', 'darkgreen'], line: {color: 'white', width: 1.5}}";
        }

        json_str += "}];\n";
        json_str += "var layout = { title: '" + locSnp.genoStr3 + "', ";
        json_str += "xaxis:{tickmode: 'array', title:'haplotype', automargin: true},";
        json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
        json_str += "barmode: 'stack'};\n";
        json_str += "Plotly.newPlot('plot_h" + divName + "', data, layout);\n";
        ofs << json_str;
    }
    //for error rate line char

    if (!locSnp.baseErrorMap.empty() && !mOptions->noErrorPlot) {
        std::vector<int> keyV;
        keyV.reserve(locSnp.baseErrorMap.size());
        std::vector<double> valueV;
        valueV.reserve(locSnp.baseErrorMap.size());
        for (const auto & posIt : locSnp.baseErrorMap) {
            keyV.push_back(posIt.first);
            valueV.push_back(posIt.second);
        }

        json_str.clear();
        json_str = "var data=[{";
        json_str += "x:[" + Stats::list2string2(keyV, keyV.size()) + "],";
        json_str += "y:[" + Stats::list2string2(valueV, valueV.size()) + "],";
        json_str += "type:'scatter', mode: 'markers', marker:{color: 'red'}, textposition: 'auto'";
        json_str += "}];\n";
        json_str += "var layout = {xaxis:{title:'Position', automargin: true},";
        json_str += "yaxis:{title:'Error rate (%)', automargin: true},";
        json_str += "shapes: [";
        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string(locSnp.aveErrorRate) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string(locSnp.aveErrorRate) + ", ";
        json_str += "line:{color: 'gray', width: 2, dash: 'dash'}";
        json_str += "}],";
        json_str += "};\n";
        json_str += "Plotly.newPlot('plot_e" + divName + "', data, layout);\n";
        ofs << json_str;
    }

    ofs << "</script>" << endl;
}

void HtmlReporter::reportSnpTablePlot(ofstream& ofs, std::string & divName, LocSnp2 & locSnp) {

    std::vector<std::string> x_vec(locSnp.ssnpsMap.size(), ""); //A
    std::vector<int> y_vec(locSnp.ssnpsMap.size(), 0);
    std::vector<double> bar_width_vec(locSnp.ssnpsMap.size(), 0.5);

    std::vector<std::string> x_vec_2(locSnp.ssnpsMap.size(), ""); //C
    std::vector<int> y_vec_2(locSnp.ssnpsMap.size(), 0);
    std::vector<double> bar_width_vec_2(locSnp.ssnpsMap.size(), 0.5);

    std::vector<std::string> x_vec_3(locSnp.ssnpsMap.size(), ""); //G
    std::vector<int> y_vec_3(locSnp.ssnpsMap.size(), 0);
    std::vector<double> bar_width_vec_3(locSnp.ssnpsMap.size(), 0.5);

    std::vector<std::string> x_vec_4(locSnp.ssnpsMap.size(), ""); //T
    std::vector<int> y_vec_4(locSnp.ssnpsMap.size(), 0);
    std::vector<double> bar_width_vec_4(locSnp.ssnpsMap.size(), 0.5);

    int i = 0;
    for (auto & it2 : locSnp.ssnpsMap) {
        x_vec[i] = x_vec_2[i] = x_vec_3[i] = x_vec_4[i] = (std::to_string(it2.first + locSnp.ft.length()) + it2.second.snp1 + "|" + it2.second.snp2);
        if (it2.second.snp1 == it2.second.snp2) {
            if (it2.second.snp1 == 'A') {
                y_vec[i] = locSnp.totHaploReads;
            } else if (it2.second.snp1 == 'C') {
                y_vec_2[i] = locSnp.totHaploReads;
            } else if (it2.second.snp1 == 'G') {
                y_vec_3[i] = locSnp.totHaploReads;
            } else if (it2.second.snp1 == 'T') {
                y_vec_4[i] = locSnp.totHaploReads;
            }
        } else {
            if (it2.second.snp1 == 'A') {
                y_vec.at(i) = locSnp.getHaploReads();

                if (it2.second.snp2 == 'C') {
                    y_vec_2[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'G') {
                    y_vec_3[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'T') {
                    y_vec_4[i] = locSnp.getHaploReads(true);
                }

            } else if (it2.second.snp1 == 'C') {
                y_vec_2.at(i) = locSnp.getHaploReads();

                if (it2.second.snp2 == 'A') {
                    y_vec[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'G') {
                    y_vec_3[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'T') {
                    y_vec_4[i] = locSnp.getHaploReads(true);
                }


            } else if (it2.second.snp1 == 'G') {
                y_vec_3.at(i) = locSnp.getHaploReads();

                if (it2.second.snp2 == 'A') {
                    y_vec[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'C') {
                    y_vec_2[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'T') {
                    y_vec_4[i] = locSnp.getHaploReads(true);
                }

            } else if (it2.second.snp1 == 'T') {
                y_vec_4.at(i) = locSnp.getHaploReads();

                if (it2.second.snp2 == 'A') {
                    y_vec[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'C') {
                    y_vec_2[i] = locSnp.getHaploReads(true);
                } else if (it2.second.snp2 == 'G') {
                    y_vec_3[i] = locSnp.getHaploReads(true);
                }
            }
        }
        i++;
    }

    ofs << "<div class='sub_section_tips'>Reads mean and thresholds for heter are in white and purple while threshold for homo loci is in yellow .</div>\n";
    ofs << "<div class='sub_section_tips'><font color='red'>Caution:</font> Position starts with <font color='red'>0</font>!</div>\n";

    ofs << "<div class='figurefull' id='plot_" + divName + "'></div>\n";

    ofs << "<div class='sub_section_tips'><font color='red'> Heter target loci are in red</font>, <font color='green'> homo target loci are in green</font>, <font color='orange'> while new heter (including homo against reference) loci are in orange</font></div>\n";
    ofs << "<div class='sub_section_tips'><span style='background-color: purple;'><font color='white'>Caution: if there are at least 2 SNPs (heter), the SNPs are still regarded as true SNPs even if the Reads ratio is between purple and yellow!</font></span></div>\n";
    ofs << "<pre overflow: scroll>\n";
    ofs << "<table class='summary_table' style='width:40%'>\n";
    ofs << "<tr style='background:#cccccc'> <td>ID</td><td>Position</td><td>Genotype</td><td>Putative Genotype</td><td>N. of reads1</td><td>N. of reads2</td><td>Reads ratio</td><td>Total reads</td></tr>\n";

    outputRow(ofs, locSnp, false);

    ofs << "</table>\n";
    ofs << "</pre>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;

    std::vector<std::string> Avec(x_vec.size(), "A");
    std::vector<std::string> Cvec(x_vec.size(), "C");
    std::vector<std::string> Gvec(x_vec.size(), "G");
    std::vector<std::string> Tvec(x_vec.size(), "T");
    string json_str = "var data=[";
    json_str += "A = {";
    json_str += "x:[" + Stats::list2string(x_vec, x_vec.size()) + "],";
    json_str += "y:[" + Stats::list2string2(y_vec, y_vec.size()) + "],";
    json_str += "text: [" + Stats::list2string(Avec, x_vec.size()) + "],";
    //json_str += "width: [" + Stats::list2string2(bar_width_vec, bar_width_vec.size()) + "],";
    json_str += "name: 'A',";
    json_str += "marker:{color:'darkgreen'},";
    json_str += "type:'bar', textposition: 'auto'";
    json_str += "},\n";

    json_str += "C = {";
    json_str += "x:[" + Stats::list2string(x_vec_2, x_vec_2.size()) + "],";
    json_str += "y:[" + Stats::list2string2(y_vec_2, y_vec_2.size()) + "],";
    json_str += "text: [" + Stats::list2string(Cvec, x_vec_2.size()) + "],";
    //json_str += "width: [" + Stats::list2string2(bar_width_vec_2, bar_width_vec_2.size()) + "],";
    json_str += "name: 'C',";
    json_str += "marker:{color:'darkred'},";
    json_str += "type:'bar', textposition: 'auto'";
    json_str += "}, \n";

    json_str += "G = {";
    json_str += "x:[" + Stats::list2string(x_vec_3, x_vec_3.size()) + "],";
    json_str += "y:[" + Stats::list2string2(y_vec_3, y_vec_3.size()) + "],";
    json_str += "text: [" + Stats::list2string(Gvec, x_vec_3.size()) + "],";
    //json_str += "width: [" + Stats::list2string2(bar_width_vec_3, bar_width_vec_3.size()) + "],";
    json_str += "name: 'G',";
    json_str += "marker:{color:'grey'},";
    json_str += "type:'bar', textposition: 'auto'";
    json_str += "}, \n";

    json_str += "T = {";
    json_str += "x:[" + Stats::list2string(x_vec_4, x_vec_4.size()) + "],";
    json_str += "y:[" + Stats::list2string2(y_vec_4, y_vec_4.size()) + "],";
    json_str += "text: [" + Stats::list2string(Tvec, x_vec_4.size()) + "],";
    //json_str += "width: [" + Stats::list2string2(bar_width_vec_4, bar_width_vec_4.size()) + "],";
    json_str += "name: 'T',";
    json_str += "marker:{color:'darkblue'},";
    json_str += "type:'bar', textposition: 'auto'";
    json_str += "}";

    json_str += "];\n";
    json_str += "var layout={autosize: true, title:'" + locSnp.name + "',";

    json_str += "shapes: [";
    json_str += "{ type: 'line', xref: 'paper', ";
    json_str += "x0: 0, ";
    json_str += "y0: " + std::to_string((double) locSnp.totHaploReads / 2) + ", ";
    json_str += "x1: 1, ";
    json_str += "y1: " + std::to_string((double) locSnp.totHaploReads / 2) + ", ";
    json_str += "line:{color: 'white', width: 4, dash: 'line'}";
    json_str += "},\n";

    if (locSnp.getNumSnps() > 1) {

        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string((double) locSnp.totHaploReads * (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH)) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string((double) locSnp.totHaploReads * (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH)) + ", ";
        json_str += "line:{color: 'goldenrod', width: 4, dash: 'line'}";
        json_str += "},\n";

        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string((double) locSnp.totHaploReads * (1 - (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH))) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string((double) locSnp.totHaploReads * (1 - (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH))) + ", ";
        json_str += "line:{color: 'goldenrod', width: 4, dash: 'line'}";
        json_str += "}\n";

    } else {

        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string((double) locSnp.totHaploReads / 2 - ((double) locSnp.totHaploReads * mOptions->mLocSnps.mLocSnpOptions.htJetter)) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string((double) locSnp.totHaploReads / 2 - ((double) locSnp.totHaploReads * mOptions->mLocSnps.mLocSnpOptions.htJetter)) + ", ";
        json_str += "line:{color: 'magenta', width: 2, dash: 'line'}";
        json_str += "},\n";

        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string((double) locSnp.totHaploReads / 2 + ((double) locSnp.totHaploReads * mOptions->mLocSnps.mLocSnpOptions.htJetter)) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string((double) locSnp.totHaploReads / 2 + ((double) locSnp.totHaploReads * mOptions->mLocSnps.mLocSnpOptions.htJetter)) + ", ";
        json_str += "line:{color: 'magenta', width: 2, dash: 'line'}";
        json_str += "},\n";

        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string((double) locSnp.totHaploReads * (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH)) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string((double) locSnp.totHaploReads * (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH)) + ", ";
        json_str += "line:{color: 'goldenrod', width: 4, dash: 'line'}";
        json_str += "},\n";

        json_str += "{ type: 'line', xref: 'paper', ";
        json_str += "x0: 0, ";
        json_str += "y0: " + std::to_string((double) locSnp.totHaploReads * (1 - (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH))) + ", ";
        json_str += "x1: 1, ";
        json_str += "y1: " + std::to_string((double) locSnp.totHaploReads * (1 - (locSnp.genoStr3 == "homo" ? mOptions->mLocSnps.mLocSnpOptions.hmPerL : mOptions->mLocSnps.mLocSnpOptions.hmPerH))) + ", ";
        json_str += "line:{color: 'goldenrod', width: 4, dash: 'line'}";
        json_str += "}\n";

    }

    json_str += "],\n";

    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(x_vec, x_vec.size()) + "],  title:'" + "SNP" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "var config = {responsive: true};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout, config);\n";
    ofs << json_str;
    ofs << "</script>" << endl;

}

void HtmlReporter::reportAllGenotype(ofstream& ofs, std::vector<std::map<std::string, std::vector<std::pair<std::string, Genotype>>>> &sortedAllGenotypeMapVec) {
    std::map<std::string, std::vector < UnitedGenotype>> uGenoMap;
    for (const auto & it : sortedAllGenotypeMapVec.at(0)) {
        std::multimap<int, Genotype> tmpGenoMMap;
        for (const auto & it2 : it.second) {
            tmpGenoMMap.insert(std::make_pair(it2.second.baseLocVar.effectiveLen, it2.second));
        }
        std::vector<UnitedGenotype> tmpUGenoVec;
        for (auto it2 = tmpGenoMMap.begin(); it2 != tmpGenoMMap.end(); it2 = tmpGenoMMap.upper_bound(it2->first)) {
            auto genoSize = it2->first;
            auto geno = tmpGenoMMap.equal_range(genoSize);
            UnitedGenotype uGeno;
            uGeno.marker = it.first;
            uGeno.effectiveLen = genoSize;

            for (std::multimap<int, Genotype>::iterator it3 = geno.first; it3 != geno.second; it3++) {
                uGeno.locVec.push_back(it3->second.baseLocVar);
                uGeno.mraSVec.push_back(it3->second.baseLocVar.mra.mStr.length());
                uGeno.numReadsVec.push_back(it3->second.numReads);
                uGeno.effectiveSeqVec.push_back(it3->second.baseLocVar.effectiveSeq.mStr);
                uGeno.seqNameVec.push_back(it3->second.baseLocVar.ff.mStr + it3->second.baseLocVar.mraName + it3->second.baseLocVar.rf.mStr);

            }
            tmpUGenoVec.push_back(uGeno);
        }
        uGenoMap[it.first] = tmpUGenoVec;
    }

    std::map<std::string, std::vector < UnitedGenotype>> uGenoMraMap;
    for (const auto & it : sortedAllGenotypeMapVec.at(1)) {
        std::multimap<int, Genotype> tmpGenoMMap;
        for (const auto & it2 : it.second) {
            tmpGenoMMap.insert(std::make_pair(it2.second.baseLocVar.mra.mStr.length(), it2.second));
        }

        std::vector<UnitedGenotype> tmpUGenoVec;
        for (auto it2 = tmpGenoMMap.begin(); it2 != tmpGenoMMap.end(); it2 = tmpGenoMMap.upper_bound(it2->first)) {
            auto mraSize = it2->first;
            auto geno = tmpGenoMMap.equal_range(mraSize);

            UnitedGenotype uGeno;
            uGeno.marker = it.first;
            uGeno.repuit = Sequence(it2->second.baseLocVar.repuit.mStr);
            uGeno.mraName = it2->second.baseLocVar.mraName;
            uGeno.mra = Sequence(it2->second.baseLocVar.mra.mStr);
            for (std::multimap<int, Genotype>::iterator it3 = geno.first; it3 != geno.second; it3++) {
                uGeno.locVec.push_back(it3->second.baseLocVar);
                uGeno.mraSVec.push_back(it3->second.baseLocVar.mra.mStr.length());
                uGeno.numReadsVec.push_back(it3->second.numReads);
            }
            tmpUGenoVec.push_back(uGeno);
        }
        uGenoMraMap[it.first] = tmpUGenoVec;
    }

    for (const auto & it : uGenoMap) {
        auto itt = sortedAllGenotypeMapVec.at(0).find(it.first);
        std::vector<int> x_vec;
        x_vec.reserve(it.second.size());
        std::vector<double> bar_width_vec;
        bar_width_vec.reserve(it.second.size());
        int nStacks = 0;
        for (const auto & it2 : it.second) {
            x_vec.push_back(it2.effectiveLen);
            if (it2.effectiveSeqVec.size() > nStacks) {
                nStacks = it2.effectiveSeqVec.size();
            }
            bar_width_vec.push_back(0.5);
        }
        std::map< std::string, std::vector<int>> stackMap;
        std::map< std::string, std::vector < std::string>> stackYlabMap;
        for (int i = 0; i < nStacks; i++) {
            std::vector<int> numReadsVec;
            numReadsVec.reserve(nStacks);
            std::vector<std::string> yLabVec;
            yLabVec.reserve(nStacks);
            for (const auto & it2 : it.second) {
                int numReads = 0;
                std::string read = "";
                if (i < it2.numReadsVec.size()) {
                    numReads = it2.numReadsVec.at(i);
                    //read = it2.seqName;
                    read = it2.seqNameVec.at(i);
                }
                numReadsVec.push_back(numReads);
                yLabVec.push_back(read);
            }
            stackMap["Genotype" + std::to_string(i)] = numReadsVec;
            stackYlabMap["Genotype" + std::to_string(i)] = yLabVec;
        }

        auto itm = sortedAllGenotypeMapVec.at(1).find(it.first);
        auto itk = uGenoMraMap.find(it.first);

        std::vector<int> xmra_vec;
        xmra_vec.reserve(itk->second.size());
        //        std::vector<std::string> ymra_label_vec;
        //        ymra_label_vec.reserve(itk->second.size());
        std::vector<double> barmra_width_vec;
        barmra_width_vec.reserve(itk->second.size());
        int nmStacks = 0;
        for (const auto & it2 : itk->second) {
            xmra_vec.push_back(it2.mra.mStr.length());
            if (it2.locVec.size() > nmStacks) {
                nmStacks = it2.locVec.size();
            }
            //ymra_label_vec.push_back(it2.mraName);
            barmra_width_vec.push_back(0.5);
        }
        std::map< std::string, std::vector<int>> stackMraMap;
        std::map< std::string, std::vector < std::string>> stackYLabMMap;
        for (int i = 0; i < nmStacks; i++) {
            std::vector<int> numReadsVec;
            numReadsVec.reserve(nStacks);
            std::vector<std::string> yLabVec;
            yLabVec.reserve(nStacks);
            for (const auto & it2 : itk->second) {
                int numReads = i < it2.numReadsVec.size() ? it2.numReadsVec.at(i) : 0;
                numReadsVec.push_back(numReads);
                std::string mra = i < it2.numReadsVec.size() ? it2.locVec.at(i).mraName : "";
                yLabVec.push_back(mra);
            }
            stackMraMap["MRA" + std::to_string(i)] = numReadsVec;
            stackYLabMMap["MRA" + std::to_string(i)] = yLabVec;
        }
        reportEachGenotype(ofs, it.first, x_vec, stackMap, stackYlabMap, bar_width_vec, itt->second,
                xmra_vec, stackMraMap, stackYLabMMap, barmra_width_vec, itm->second);
        x_vec.clear();
        x_vec.shrink_to_fit();
        bar_width_vec.clear();
        bar_width_vec.shrink_to_fit();
        stackMap.clear();
        stackYlabMap.clear();

        xmra_vec.clear();
        xmra_vec.shrink_to_fit();
        barmra_width_vec.clear();
        barmra_width_vec.shrink_to_fit();
        stackMraMap.clear();
        stackYLabMMap.clear();
    }
    uGenoMap.clear();
    uGenoMraMap.clear();
    sortedAllGenotypeMapVec.clear();
}

void HtmlReporter::reportEachGenotype(ofstream& ofs, std::string marker,
        std::vector<int> & x_vec, std::map< std::string, std::vector<int>> &stackMap,
        std::map< std::string, std::vector<std::string>> &stackYlabMap, std::vector<double> & bar_width_vec,
        std::vector<std::pair<std::string, Genotype>> &outGenotype,
        std::vector<int> & xmra_vec, std::map< std::string, std::vector<int>> &stackMraMap,
        std::map< std::string, std::vector<std::string>> &stackYLabMMap,
        std::vector<double> & barmra_width_vec,
        std::vector<std::pair<std::string, Genotype>> &outGenotypeMra) {
    std::string subsection = "Marker: " + marker;
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title = marker;

    ofs << "<div class='subsection_title' onclick=showOrHide('" + divName + "')><a name='" + subsection + "'>" + subsection + "<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    //ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" + divName + "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each allele size will be shown on mouse over.</div>\n";

    ofs << "<div class='left'>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    //    ofs << "<div class='right'>\n";
    //    ofs << "<div class='figure' id='plot_mra" + divName + "'></div>\n";
    //    ofs << "</div>";
    //ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";

    ofs << "<div class='sub_section_tips'>SNPs/artifacts are highlighted in red. N. of Reads are in red are the warnings</div>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr style='background:#cccccc'> <td>Marker</td><td>Repeat unit</td><td>MRA base</td><td>MRA size</td><td>Allele size</td><td>N. of Reads</td><td align = 'right'>Forward flanking region</td><td align='center'>MRA</td><td align='left'>Reverse flanking region</td></tr>\n";
    auto it = mOptions->mLocVars.refLocMap.find(marker);
    if (it != mOptions->mLocVars.refLocMap.end()) {
        ofs << "<tr style='color:blue'>";
        ofs << "<td>Reference</td>" <<
                "<td>" + it->second.repuitAll.mStr + "</td>" <<
                "<td>" + it->second.mraBase + "</td>" <<
                "<td>" + std::to_string(it->second.mra.length()) + "</td>" <<
                "<td>" + std::to_string(it->second.effectiveLen) + "</td>" <<
                "<td>N.A.</td>" <<
                "<td align='right'>" + highligher(it->second.ff.mStr, it->second.refSnpsSetffMap[basename(mOptions->prefix)]) + "</td>" << //could use <xmp>
                "<td align='center'>" + it->second.mraName + "</td>" <<
                "<td align='left'>" + highligher(it->second.rf.mStr, it->second.refSnpsSetrfMap[basename(mOptions->prefix)]) + "</td>";
        ofs << "</tr>";
        outputRow(ofs, marker, outGenotypeMra, mOptions);
    }

    ofs << "</table>\n";

    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;

    string json_str = "";
    for (auto & it : stackMap) {
        auto itt = stackYlabMap.find(it.first);
        json_str += "var " + it.first + " = {";
        json_str += "x:[" + Stats::list2string(x_vec, x_vec.size()) + "],";
        json_str += "y:[" + Stats::list2string2(it.second, it.second.size()) + "],";
        //json_str += "text: [" + Stats::list2string(itt->second, itt->second.size()) + "],";
        //json_str += "width: [" + Stats::list2string2(barmra_width_vec, barmra_width_vec.size()) + "],";
        json_str += "type:'bar', textposition: 'auto',";
        json_str += "};\n";
    }

    json_str += "var data = [" + Stats::list2string(stackMap, stackMap.size()) + "];\n";
    json_str += "var layout={title:'" + title + "',";
    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(x_vec, x_vec.size()) + "],  title:'" + "Allele size (bp)" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    for (auto & it : stackMraMap) {
        auto itt = stackYLabMMap.find(it.first);
        json_str += "var " + it.first + " = {";
        json_str += "x:[" + Stats::list2string(xmra_vec, xmra_vec.size()) + "],";
        json_str += "y:[" + Stats::list2string2(it.second, it.second.size()) + "],";
        //json_str += "text: [" + Stats::list2string(itt->second, itt->second.size()) + "],";
        //json_str += "width: [" + Stats::list2string(barmra_width_vec, barmra_width_vec.size()) + "],";
        json_str += "type:'bar', textposition: 'auto',";
        json_str += "};\n";
    }

    json_str += "var data = [" + Stats::list2string(stackMraMap, stackMraMap.size()) + "];\n";
    json_str += "var layout={title:'" + title + "',";
    json_str += "xaxis:{tickmode: 'array', tickvals:[" + Stats::list2string(xmra_vec, xmra_vec.size()) + "],  title:'" + "MRA size (bp)" + "', automargin: true},";
    json_str += "yaxis:{title:'Number of reads', automargin: true}, ";
    json_str += "barmode: 'stack'};\n";
    json_str += "Plotly.newPlot('plot_mra" + divName + "', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
}

std::string HtmlReporter::highligher(std::string & str, std::map<int, std::string> & snpsMap) {
    if (snpsMap.empty()) {
        return str;
    }
    std::string hstr = str;
    for (auto it = snpsMap.rbegin(); it != snpsMap.rend(); it++) {
        hstr.insert(it->first + 1, "</mark>");
        hstr.insert(it->first, "<mark>");
    }
    return hstr;
}

std::string HtmlReporter::highligher(std::string & str, std::set<int> & snpsSet) {
    if (snpsSet.empty()) {
        return str;
    }
    std::string hstr = str;
    for (auto it = snpsSet.rbegin(); it != snpsSet.rend(); it++) {
        hstr.insert(*it + 1, "</mark>");
        hstr.insert(*it, "<mark>");
    }
    return hstr;
}

std::string HtmlReporter::highligher(LocSnp2 & locSnp, bool ref, std::string refStr, std::string tarStr, std::set<int> & posSet) {
    if (ref) {
        std::string ft = locSnp.ft.mStr;
        for (int i = ft.length() - 1; i >= 0; i--) {
            ft.insert(i + 1, "</mark4>");
            ft.insert(i, "<mark4>");
        }

        std::string rt = locSnp.rt.mStr;
        for (int i = rt.length() - 1; i >= 0; i--) {
            rt.insert(i + 1, "</mark4>");
            rt.insert(i, "<mark4>");
        }

        std::string hstr = locSnp.ref.mStr;
        for (auto it = locSnp.snpPosSet.rbegin(); it != locSnp.snpPosSet.rend(); it++) {
             if (locSnp.snpPosSetHaplo.find(*it) != locSnp.snpPosSetHaplo.end()) {
                if (locSnp.refSnpPosSet.find(*it) != locSnp.refSnpPosSet.end()) {
                    if(locSnp.ref.mStr[*it] == locSnp.ssnpsMap[*it].snp1 && locSnp.ref.mStr[*it] == locSnp.ssnpsMap[*it].snp2) {
                        hstr.insert(*it + 1, "</mark3>");
                        hstr.insert(*it, "<mark3>");
                    } else {
                        hstr.insert(*it + 1, "</mark>");
                        hstr.insert(*it, "<mark>");
                    }
                } else {
                    hstr.insert(*it + 1, "</mark2>");
                    hstr.insert(*it, "<mark2>");
                }
            } else {
                hstr.insert(*it + 1, "</mark4>");
                hstr.insert(*it, "<mark4>");
            }
        }
        return (ft + hstr + rt);

    } else {

        std::string ft(locSnp.ft.length(), '*');
        std::string rt(locSnp.rt.length(), '*');

        if (posSet.empty()) {
            return (ft + tarStr + rt);
        } else {
            std::string hstr = tarStr;
            //            for (int i = (tarStr.length() -1); i >= 0; i--) {
            //                if (posSet.find(i) == posSet.end()) {
            //                    
            //                    switch(tarStr[i]){
            //                        case 'A':
            //                            hstr.insert(i + 1, "</fontA>");
            //                            hstr.insert(i, "<fontA>");
            //                            break;
            //                        case 'C':
            //                            hstr.insert(i + 1, "</fontC>");
            //                            hstr.insert(i, "<fontC>");
            //                            break;
            //                        case 'G':
            //                            hstr.insert(i + 1, "</fontG>");
            //                            hstr.insert(i, "<fontG>");
            //                            break;
            //                        case 'T':
            //                            hstr.insert(i + 1, "</fontT>");
            //                            hstr.insert(i, "<fontT>");
            //                            break;
            //                        default:
            //                            break;
            //                    }
            //
            //                } else {
            //                    if (refStr[i] != tarStr[i]) {
            //                        if (locSnp.refSnpPosSet.find(i) == locSnp.refSnpPosSet.end()) {
            //                            hstr.insert(i + 1, "</mark2>");
            //                            hstr.insert(i, "<mark2>");
            //                        } else {
            //                            hstr.insert(i + 1, "</mark>");
            //                            hstr.insert(i, "<mark>");
            //                        }
            //                    } else {
            //                        switch (tarStr[i]) {
            //                            case 'A':
            //                                hstr.insert(i + 1, "</fontA>");
            //                                hstr.insert(i, "<fontA>");
            //                                break;
            //                            case 'C':
            //                                hstr.insert(i + 1, "</fontC>");
            //                                hstr.insert(i, "<fontC>");
            //                                break;
            //                            case 'G':
            //                                hstr.insert(i + 1, "</fontG>");
            //                                hstr.insert(i, "<fontG>");
            //                                break;
            //                            case 'T':
            //                                hstr.insert(i + 1, "</fontT>");
            //                                hstr.insert(i, "<fontT>");
            //                                break;
            //                            default:
            //                                break;
            //                        }
            //                    }
            //                }
            //            }

            for (auto it = posSet.rbegin(); it != posSet.rend(); it++) {
                if (locSnp.snpPosSetHaplo.find(*it) != locSnp.snpPosSetHaplo.end()) {
                    if (locSnp.refSnpPosSet.find(*it) != locSnp.refSnpPosSet.end()) {
                        if(refStr[*it] == tarStr[*it]) {
                            hstr.insert(*it + 1, "</mark3>");
                            hstr.insert(*it, "<mark3>");
                        } else {
                            hstr.insert(*it + 1, "</mark>");
                            hstr.insert(*it, "<mark>");
                        }
                    } else {
                        hstr.insert(*it + 1, "</mark2>");
                        hstr.insert(*it, "<mark2>");
                    }
                } else {
                    hstr.insert(*it + 1, "</mark4>");
                    hstr.insert(*it, "<mark4>");
                }
            }
            return (ft + hstr + rt);
        }
    }
}

void HtmlReporter::printHeader(ofstream& ofs) {
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>Seq2Sat report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printCSS(ofstream& ofs) {
    ofs << "<style type=\"text/css\">" << endl;
    //ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "td {border:1px; solid #dddddd;padding:5px;font-size:12px; width:1px; white-space:nowrap; border:1px solid gray;}" << endl;
    //ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2px;border-collapse:collapse; table-layout:auto; border:1px solid gray;}" << endl;
    ofs << ".col1 {width:320px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:auto; height:auto;}" << endl;
    ofs << ".figurefull {width:80%; height:auto;}" << endl;
    ofs << ".figurehalf {width:50%; height:auto;}" << endl;
    //ofs << ".figure {width:800px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".sub_section_div {font-size:13px;padding-left:10px;text-align:left; margin-top:10px;}" << endl;
    ofs << ".sub_section_title {color:#ffffff;font-size:13px;padding-left:10px;text-align:left;background:#009900; margin-top:10px; width:20%;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:14px;padding:7px;text-align:left;background:#009900; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:12px;padding:12px;margin-top:10px;text-align:left;color:blue}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#009900;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:12px;padding-bottom:3px;text-align:left;}" << endl;
    ofs << ".left, .right{display: inline-block}" << endl;
    ofs << "mark{background-color: red; color: white;}" << endl;
    ofs << "mark2{background-color: orange; color: white;}" << endl;
    ofs << "mark3{background-color: green; color: white;}" << endl;
    ofs << "mark4{background-color: gray; color: white;}" << endl;
    ofs << "fontA{color: green;}" << endl;
    ofs << "fontC{color: red;}" << endl;
    ofs << "fontG{color: gray;}" << endl;
    ofs << "fontT{color: blue;}" << endl;
    ofs << "pre{overflow: auto; width:0; min-width:100%;}" << endl;
    ofs << "</style>" << endl;
}

void HtmlReporter::printJS(ofstream& ofs) {
    ofs << "<script src='https://cdn.plot.ly/plotly-2.23.2.min.js' charset='utf-8'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }\n" << endl;
    ofs << "</script>" << endl;
}

const string HtmlReporter::getCurrentSystemTime() {
    auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm* ptm = localtime(&tt);
    char date[60] = {0};
    sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
            (int) ptm->tm_year + 1900, (int) ptm->tm_mon + 1, (int) ptm->tm_mday,
            (int) ptm->tm_hour, (int) ptm->tm_min, (int) ptm->tm_sec);
    return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs) {
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>" << command << "</p>";
    ofs << "Seq2Sat " << SEQ2SAT_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}
