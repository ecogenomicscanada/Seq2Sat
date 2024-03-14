#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <functional>
#include "fastareader.h"

Options::Options(){
    in1 = "";
    in2 = "";
    out1 = "";
    out2 = "";
    outFR = false;
    outFRFile = "";
    reportTitle = "seq2sat report";
    thread = 1;
    compression = 2;
    phred64 = false;
    dontOverwrite = false;
    inputFromSTDIN = false;
    outputToSTDOUT = false;
    readsToProcess = 0;
    interleavedInput = false;
    insertSizeMax = 512;
    overlapRequire = 30;
    overlapDiffLimit = 5;
    overlapDiffPercentLimit = 20;
    verbose = false;
    seqLen1 = 151;
    seqLen2 = 151;
    mergerOverlappedPE = true;
    prefix = "";
    debug = false;
    var = "";
    locFile = "";
    sexFile = "";
    sampleTable = "";
    samples.clear();
    minAmpliconEffectiveLen = 6;
    revCom = false;
    noErrorPlot = false;
    noSnpPlot = false;
    noPlot = false;
}

void Options::init() {
}

bool Options::isPaired() {
    return in2.length() > 0 || interleavedInput;
}

bool Options::adapterCuttingEnabled() {
    if(adapter.enabled){
        if(isPaired() || !adapter.sequence.empty())
            return true;
    }
    return false;
}

bool Options::polyXTrimmingEnabled() {
    return polyXTrim.enabled;
}

void Options::loadFastaAdapters() {
    if(adapter.fastaFile.empty()) {
        adapter.hasFasta = false;
        return;
    }

    check_file_valid(adapter.fastaFile);

    FastaReader reader(adapter.fastaFile);
    reader.readAll();

    map<string, string> contigs = reader.contigs();
    map<string, string>::iterator iter;
    for(iter = contigs.begin(); iter != contigs.end(); iter++) {
        if(iter->second.length()>=6) {
            adapter.seqsInFasta.push_back(iter->second);
        }
        else {
            cerr << "skip too short adapter sequence in " <<  adapter.fastaFile << " (6bp required): " << iter->second << endl;
        }
    }

    if(adapter.seqsInFasta.size() > 0) {
        adapter.hasFasta = true;
    } else {
        adapter.hasFasta = false;
    }
}

bool Options::validate() {
    if(in1.empty()) {
        if(!in2.empty())
            error_exit("read2 input is specified by <in2>, but read1 input is not specified by <in1>");
        if(inputFromSTDIN)
            in1 = "/dev/stdin";
        else
            error_exit("read1 input should be specified by --in1, or enable --stdin if you want to read STDIN");
    } else {
        check_file_valid(in1);
    }
    
    if(prefix.empty()){
        prefix = getPrefix(in1);
    }
    
    htmlFile = prefix + ".html";
    jsonFile = prefix + ".json";

    if(!in2.empty()) {
        check_file_valid(in2);
    }

    // if output to STDOUT, then...
    if(outputToSTDOUT) {
        if(isPaired())
            cerr << "Enable interleaved output mode for paired-end input." << endl;
    }

    if(in2.empty() && !interleavedInput && !out2.empty()) {
        error_exit("read2 output is specified (--out2), but neighter read2 input is not specified (--in2), nor read1 is interleaved.");
    }

    if(!in2.empty() || interleavedInput) {
        if(!out1.empty() && out2.empty()) {
            error_exit("paired-end input, read1 output should be specified together with read2 output (--out2 needed) ");
        }
        if(out1.empty() && !out2.empty()) {
            error_exit("paired-end input, read1 output should be specified (--out1 needed) together with read2 output ");
        }
    }

    if(!in2.empty() && interleavedInput) {
        error_exit("<in2> is not allowed when <in1> is specified as interleaved mode by (--interleaved_in)");
    }

    if(!out1.empty()) {
        //check_file_writable(out1);
        if(out1 == out2) {
            error_exit("read1 output (--out1) and read1 output (--out2) should be different");
        }
        if(dontOverwrite && file_exists(out1)) {
            error_exit(out1 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    if(!out2.empty()) {
        //check_file_writable(out2);
        if(dontOverwrite && file_exists(out2)) {
            error_exit(out2 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(dontOverwrite) {
        if(file_exists(jsonFile)) {
            error_exit(jsonFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(file_exists(htmlFile)) {
            error_exit(htmlFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(outFR){
        if(outFRFile.empty()){
            error_exit(prefix + " is empty; please specify the prefix");
        }
    }
    
    if(compression < 1 || compression > 9)
        error_exit("compression level (--compression) should be between 1 ~ 9, 1 for fastest, 9 for smallest");

    if(readsToProcess < 0)
        error_exit("the number of reads to process (--reads_to_process) cannot be negative");

    if(thread < 1) {
        thread = 1;
    } else if(thread > 16) {
        cerr << "WARNING: seq2sat uses up to 16 threads although you specified " << thread << endl;
        thread = 16;
    }

    if(trim.front1 < 0 || trim.front1 > 30)
        error_exit("trim_front1 (--trim_front1) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail1 < 0 || trim.tail1 > 100)
        error_exit("trim_tail1 (--trim_tail1) should be 0 ~ 100, suggest 0 ~ 4");

    if(trim.front2 < 0 || trim.front2 > 30)
        error_exit("trim_front2 (--trim_front2) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail2 < 0 || trim.tail2 > 100)
        error_exit("trim_tail2 (--trim_tail2) should be 0 ~ 100, suggest 0 ~ 4");

    if(qualfilter.qualifiedQual - 33 < 0 || qualfilter.qualifiedQual - 33 > 93)
        error_exit("qualitified phred (--qualified_quality_phred) should be 0 ~ 93, suggest 10 ~ 20");

    if(qualfilter.avgQualReq < 0 || qualfilter.avgQualReq  > 93)
        error_exit("average quality score requirement (--average_qual) should be 0 ~ 93, suggest 20 ~ 30");

    if(qualfilter.unqualifiedPercentLimit < 0 || qualfilter.unqualifiedPercentLimit > 100)
        error_exit("unqualified percent limit (--unqualified_percent_limit) should be 0 ~ 100, suggest 20 ~ 60");

    if(qualfilter.nBaseLimit < 0 || qualfilter.nBaseLimit > 50)
        error_exit("N base limit (--n_base_limit) should be 0 ~ 50, suggest 3 ~ 10");

    if(lengthFilter.requiredLength < 0 )
        error_exit("length requirement (--length_required) should be >0, suggest 15 ~ 100");

    if(qualityCut.enabledFront || qualityCut.enabledTail || qualityCut.enabledRight) {
        if(qualityCut.windowSizeShared < 1 || qualityCut.windowSizeShared > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_window_size) should be between 1~1000.");
        if(qualityCut.qualityShared < 1 || qualityCut.qualityShared > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeFront < 1 || qualityCut.windowSizeFront > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_front_window_size) should be between 1~1000.");
        if(qualityCut.qualityFront < 1 || qualityCut.qualityFront > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_front_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeTail < 1 || qualityCut.windowSizeTail > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_tail_window_size) should be between 1~1000.");
        if(qualityCut.qualityTail < 1 || qualityCut.qualityTail > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_tail_mean_quality) should be 1 ~ 30, suggest 13 ~ 20.");
        if(qualityCut.windowSizeRight < 1 || qualityCut.windowSizeRight > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_right_window_size) should be between 1~1000.");
        if(qualityCut.qualityRight < 1 || qualityCut.qualityRight > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_right_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
    }

    if(adapter.sequence!="auto" && !adapter.sequence.empty()) {
        // validate adapter sequence for single end adapter trimming
        if(adapter.sequence.length() <= 3)
            error_exit("the sequence of <adapter_sequence> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequence.length(); i++) {
            char c = adapter.sequence[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence> can only have bases in {A, T, C, G}, but the given sequence is: " + adapter.sequence);
            }
        }

        adapter.hasSeqR1 = true;
    }

    if(adapter.sequenceR2!="auto" && !adapter.sequenceR2.empty()) {
        // validate adapter sequenceR2 for single end adapter trimming
        if(adapter.sequenceR2.length() <= 3)
            error_exit("the sequence of <adapter_sequence_r2> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequenceR2.length(); i++) {
            char c = adapter.sequenceR2[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence_r2> can only have bases in {A, T, C, G}, but the given sequenceR2 is: " + adapter.sequenceR2);
            }
        }

        adapter.hasSeqR2 = true;
    }

    if(correction.enabled && !isPaired()) {
        cerr << "WARNING: base correction is only appliable for paired end data, ignoring -c/--correction" << endl;
        correction.enabled = false;
    }

    if(umi.enabled) {
        if(umi.location == UMI_LOC_READ1 || umi.location == UMI_LOC_READ2 || umi.location == UMI_LOC_PER_READ) {
            if(umi.length<1 || umi.length>100)
                error_exit("UMI length should be 1~100");
            if(umi.skip<0 || umi.skip>100)
                error_exit("The base number to skip after UMI <umi_skip> should be 0~100");
        }else {
            if(umi.skip>0)
                error_exit("Only if the UMI location is in read1/read2/per_read, you can skip bases after UMI");
            if(umi.length>0)
                error_exit("Only if the UMI location is in read1/read2/per_read, you can set the UMI length");
        }
        if(!umi.prefix.empty()) {
            if(umi.prefix.length() >= 10)
                error_exit("UMI prefix should be shorter than 10");
            for(int i=0; i<umi.prefix.length(); i++) {
                char c = umi.prefix[i];
                if( !(c>='A' && c<='Z') && !(c>='a' && c<='z') && !(c>='0' && c<='9')) {
                    error_exit("UMI prefix can only have characters and numbers, but the given is: " + umi.prefix);
                }
            }
        }
        if(!umi.separator.empty()) {
            if(umi.separator.length()>10)
                error_exit("UMI separator cannot be longer than 10 base pairs");
            // validate bases
            for(int i=0; i<umi.separator.length(); i++) {
                char c = umi.separator[i];
                if(c!='A' && c!='T' && c!='C' && c!='G') {
                    error_exit("UMI separator can only have bases in {A, T, C, G}, but the given sequence is: " + umi.separator);
                }
            }
        }

    }
    
    
    if(mLocSnps.mLocSnpOptions.minSeqsPer > 1 || mLocSnps.mLocSnpOptions.minSeqsPer <= 0) error_exit("--minSeqsPerSnp must be < 100 and > 0");
    if(mLocSnps.mLocSnpOptions.htJetter >= 0 
            && ((0.5 + mLocSnps.mLocSnpOptions.htJetter) < mLocSnps.mLocSnpOptions.hmPerL )
            && mLocSnps.mLocSnpOptions.hmPerL < mLocSnps.mLocSnpOptions.hmPerH 
            && mLocSnps.mLocSnpOptions.hmPerH <= 1){
    } else {
        error_exit("--minSeqsPerSnp must be < 100");
    }

    if (locFile.empty()) {
        error_exit("locus file is empty, please provide a valid file!");
    } else {
        check_file_valid(locFile);
    }

    if (var == "ssr") {
        mVarType = ssr;
    } else if (var == "snp") {
        mVarType = snp;
    } else {
        error_exit("Invalid variance type, please specify a variance type using --var");
    }

    if(mEdOptions.mode == "HW"){
        mEdOptions.modeCode = EDLIB_MODE_HW;
    } else if(mEdOptions.mode == "SHW"){
        mEdOptions.modeCode = EDLIB_MODE_SHW;
    } else if(mEdOptions.mode == "NW"){
        mEdOptions.modeCode = EDLIB_MODE_NW;
    } else {
        error_exit("Invalid mode, please specify a mode using --mode");
    }
    
    if(mEdOptions.findStartLocation) {
        mEdOptions.alignTask = EDLIB_TASK_LOC;
    } else if (mEdOptions.findAlignment){
        mEdOptions.alignTask = EDLIB_TASK_PATH;
    } else {
        mEdOptions.alignTask = EDLIB_TASK_DISTANCE;
    }
    
    if(mLocSnps.mLocSnpOptions.maxRows4Align < 2){
        cerr << "maxRows4Align " << mLocSnps.mLocSnpOptions.maxRows4Align << " is less than the minimum 2, and it is going to be set to 2!" << endl;
        mLocSnps.mLocSnpOptions.maxRows4Align = 2;
    }
    
    return true;
}

bool Options::shallDetectAdapter(bool isR2) {
    if(!adapter.enabled)
        return false;

    if(isR2) {
        return isPaired() && adapter.detectAdapterForPE && adapter.sequenceR2 == "auto";
    } else {
        if(isPaired())
            return adapter.detectAdapterForPE && adapter.sequence == "auto";
        else
            return adapter.sequence == "auto";
    }
}

string Options::getAdapter1(){
    if(adapter.sequence == "" || adapter.sequence == "auto")
        return "unspecified";
    else
        return adapter.sequence;
}

string Options::getAdapter2(){
    if(adapter.sequenceR2 == "" || adapter.sequenceR2 == "auto")
        return "unspecified";
    else
        return adapter.sequenceR2;
}

void Options::readLocFile(){

    if (locFile.empty()) {
        error_exit("locus file is empty, please provide a valid file!");
    } else {
        check_file_valid(locFile);
    }
    
    if(verbose){
        loginfo("Reading locus file: " + locFile);
    }
    
    std::ifstream fileIn;
    fileIn.open(locFile.c_str());
    if(!fileIn.is_open()){
        error_exit("Can not open locus File: " + locFile);
    };
    
    const int maxLine = 10000;
    char line[maxLine];
    int readed = 0;
    std::vector<std::string> splitVec;
    std::string lineStr;

    if (mVarType == ssr) {
        mLocVars.refLocMap.clear();
        splitVec.reserve(8);
        while (fileIn.getline(line, maxLine)) {
            readed = strlen(line);
            if (line[readed - 1] == '\n' || line[readed - 1] == '\r') {
                line[readed - 1] = '\0';
                if (line[readed - 2] == '\r') {
                    line[readed - 2] = '\0';
                }
            }
            lineStr = std::string(line);
            splitStr(lineStr, splitVec);
            if (splitVec.size() == 8) {
                LocVar tmpLocVar;
                tmpLocVar.name = splitVec[0];
                tmpLocVar.fp = splitVec[1];
                tmpLocVar.rp = revCom ? Sequence(splitVec[2]).reverseComplement().mStr : splitVec[2];
                if (splitVec[3] != "NA") {
                    tmpLocVar.ff = splitVec[3];
                } else {
                    tmpLocVar.ff = Sequence("");
                }

                if (splitVec[4] != "NA") {
                    tmpLocVar.rf = splitVec[4];
                } else {
                    tmpLocVar.rf = Sequence("");
                }
                tmpLocVar.repuitAll = splitVec[5];
                if (splitVec[5].find("|") != std::string::npos) {
                    std::vector<std::string> mrasVec;
                    splitStr(splitVec[5], mrasVec, "|");
                    tmpLocVar.repuit = mrasVec.at(0);
                    tmpLocVar.repuit2 = mrasVec.at(1);
                    if(tmpLocVar.repuit.length() == tmpLocVar.repuit2.length()){
                        tmpLocVar.repuitAllLen = tmpLocVar.repuit.length();
                    } else {
                        tmpLocVar.repuitAllLen = tmpLocVar.repuit.length() + tmpLocVar.repuit2.length();
                    }
                } else {
                    tmpLocVar.repuit = splitVec[5];
                    tmpLocVar.repuitAllLen = 0;
                }
                tmpLocVar.nRep = std::stoi(splitVec[6]);
                if (tmpLocVar.repuit.mStr.length() == 2 || tmpLocVar.repuit.mStr.length() == 3) {
                    tmpLocVar.edCutoff = 1;
                } else if (tmpLocVar.repuit.mStr.length() == 4 || tmpLocVar.repuit.mStr.length() == 5) {
                    tmpLocVar.edCutoff = 2;
                } else if (tmpLocVar.repuit.mStr.length() == 6 || tmpLocVar.repuit.mStr.length() == 7) {
                    tmpLocVar.edCutoff = 3;
                } else if (tmpLocVar.repuit.mStr.length() == 8 || tmpLocVar.repuit.mStr.length() == 9) {
                    tmpLocVar.edCutoff = 4;
                } else {
                    tmpLocVar.edCutoff = 5;
                }

                tmpLocVar.mra = splitVec[7];
                tmpLocVar.effectiveSeq = Sequence(tmpLocVar.ff.mStr + tmpLocVar.mra.mStr + tmpLocVar.rf.mStr);
                tmpLocVar.locSeq = Sequence(std::string(tmpLocVar.fp.mStr + tmpLocVar.effectiveSeq.mStr + tmpLocVar.rp.mStr));
                tmpLocVar.locLen = tmpLocVar.locSeq.length();
                tmpLocVar.effectiveLen = tmpLocVar.effectiveSeq.length();
                tmpLocVar.mraName = getGenotype(tmpLocVar.mra.mStr, tmpLocVar.repuit.mStr, tmpLocVar.repuit2.mStr);
                tmpLocVar.mraBase = getMraBase(tmpLocVar.mraName);
                //cCout(tmpLocVar.mraName, 'r');
                mLocVars.refLocMap[tmpLocVar.name] = tmpLocVar;
            }
            splitVec.clear();
        } 
    } else if(mVarType == snp) {
        mLocSnps.refLocMap.clear();
        splitVec.reserve(7);
        std::vector<std::string> posVec;
        while(fileIn.getline(line, maxLine)) {
            readed = strlen(line);
            if(readed == 0) continue;
            if (line[readed - 1] == '\n' || line[readed - 1] == '\r') {
                line[readed - 1] = '\0';
                if (line[readed - 2] == '\r') {
                    line[readed - 2] = '\0';
                }
            }
            lineStr = std::string(line);
            if(lineStr.empty()) continue;
            splitStr(lineStr, splitVec);
            if(splitVec.size() == 7){
                LocSnp2 tmpLocSnp;
                tmpLocSnp.name = splitVec[0];
                tmpLocSnp.fp = splitVec[1];
                tmpLocSnp.rp = revCom ? Sequence(splitVec[2]).reverseComplement().mStr : splitVec[2];
                tmpLocSnp.trimPos = std::make_pair(std::stoi(splitVec[3]), std::stoi(splitVec[4]));
                
                if (!(tmpLocSnp.trimPos.first >= 0 &&
                        tmpLocSnp.trimPos.first < (splitVec[6].length() - tmpLocSnp.trimPos.second) &&
                        tmpLocSnp.trimPos.second >= 0 &&
                        tmpLocSnp.trimPos.second < (splitVec[6].length() - tmpLocSnp.trimPos.first))) {
                    error_exit(tmpLocSnp.name  + " trimming site is not correct!");
                }
                
                tmpLocSnp.ft = Sequence(splitVec[6].substr(0, tmpLocSnp.trimPos.first));
                tmpLocSnp.rt = Sequence(splitVec[6].substr((splitVec[6].length() - tmpLocSnp.trimPos.second)));
                tmpLocSnp.ref = splitVec[6].substr(tmpLocSnp.trimPos.first, splitVec[6].length() - tmpLocSnp.trimPos.second - tmpLocSnp.trimPos.first);
                
                if (splitVec[5].find("|") != std::string::npos) {
                    splitStr(splitVec[5], posVec, "|");
                } else {
                    posVec.push_back(splitVec[5]);
                }

                for (auto & itt : posVec) {
                    if (itt == "NA") break;
                    auto pos = std::stoi(itt);
                    if (pos > tmpLocSnp.trimPos.first && pos < (tmpLocSnp.trimPos.first + tmpLocSnp.ref.length())) {
                        pos -= tmpLocSnp.trimPos.first;
                        tmpLocSnp.refSnpPosSet.insert(pos);
                        tmpLocSnp.snpPosSetHaplo.insert(pos);
                        //tmpLocSnp.totPosSet.insert(pos);
                    }
                }
                posVec.clear();
                mLocSnps.refLocMap[tmpLocSnp.name] = tmpLocSnp;
                //tmpLocSnp.print();
            } else {
                error_exit("Your locus " + lineStr + " does not have 7 columns!");
            }
            splitVec.clear();
        }       
    } else {
        error_exit("You have to use --var to specify ssr or snp");
    }
    
    fileIn.close();
    if (verbose) {
        loginfo("Read loci of: " + std::to_string(mVarType == ssr ? mLocVars.refLocMap.size() : mLocSnps.refLocMap.size()));
    }
}

void Options::readSexLoc(){
    if (sexFile.empty()) {
        return;
    } else {
        check_file_valid(sexFile);
    }
    
    if(verbose){
        loginfo("Reading sex file: " + sexFile);
    }
    
    std::ifstream fileIn;
    fileIn.open(sexFile.c_str());
    if(!fileIn.is_open()){
        error_exit("Can not open sex loci File: " + sexFile);
    };
    
    const int maxLine = 10000;
    char line[maxLine];
    int readed = 0;
    std::vector<std::string> splitVec;
    std::string lineStr;
    splitVec.reserve(5);
        while (fileIn.getline(line, maxLine)) {
        readed = strlen(line);
        if (line[readed - 1] == '\n' || line[readed - 1] == '\r') {
            line[readed - 1] = '\0';
            if (line[readed - 2] == '\r') {
                line[readed - 2] = '\0';
            }
        }
        lineStr = std::string(line);
        splitStr(lineStr, splitVec);
        if (splitVec.size() == 5) {
            mSex.sexMarker = splitVec[0];
            mSex.primerF = Sequence(splitVec[1]);
            mSex.primerR = Sequence(splitVec[2]);
            mSex.refX = Sequence(splitVec[3]);
            mSex.refY = Sequence(splitVec[4]);
            mSex.lengthEqual = mSex.refX.length() == mSex.refY.length() ? true : false;
        }
        splitVec.clear();
    } 
    
    fileIn.close();
}

void Options::parseSampleTable(){
    if(sampleTable.empty()) error_exit("sample table file should be specified by --sampleTable");

    check_file_valid(sampleTable);
    if (verbose) {
        std::string msg = "Reading sample table from file " + sampleTable;
        loginfo(msg);
    }
    
    ifstream infile;
    infile.open(sampleTable.c_str(), ifstream::in);
    const int maxLine = 1000;
    char line[maxLine];
    std::vector<std::string> spltVec;
    while(infile.getline(line, maxLine)){
        int readed = strlen(line);
        if(readed >= 2){
            if(line[readed - 1] == '\n' || line[readed - 1] == '\r'){
                line[readed - 1] = '\0';
                if(line[readed - 2] == '\r'){
                    line[readed - 2] = '\0';
                }
            }
        }
        
        std::string linestr(line);
        spltVec.clear();
        splitStr(linestr, spltVec, "\t");
        
        Sample s;
        s.prefix = trimStr(spltVec[0]);
        s.path = dirname(s.prefix);
        s.in1 = trimStr(spltVec[1]);
        if(spltVec.size() == 3){
            s.feature = trimStr(spltVec[2]);
        } else if(spltVec.size() == 4){
            s.in2 = trimStr(spltVec[2]);
            s.feature = trimStr(spltVec[3]);
        } else {
            error_exit("sample table must be 2 columns with sample name, forward reads or 3 columns with sample name, forward reads, reverse reads");
        }
        if(debug) cCout("reading meta: " + s.prefix, 'g');
        this->samples.push_back(s);
    }
    infile.close();
}
