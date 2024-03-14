#include <set>

#include "genotype.h"

LocVar::LocVar() {
    this->name = "";
    this->fp = Sequence("");
    this->rp = Sequence("");
    this->ff = Sequence("");
    this->rf = Sequence("");
    this->repuitAll = Sequence("");
    this->repuit = Sequence("");
    this->repuit2 = Sequence("");
    this->nRep = 0;
    this->edCutoff = 0;
    this->mra = Sequence("");
    this->locSeq = Sequence("");
    this->effectiveSeq = Sequence("");
    this->locLen = 0;
    this->effectiveLen = 0;
    this->mraName = "";
    this->snpsMapff.clear();
    this->snpsMaprf.clear();
    this->refSnpsSetffMap.clear();
    this->refSnpsSetrfMap.clear();
    this->mraBase = "";
    this->totalReads = 0;
    this->repuitAllLen = 0;
}

LocVar::LocVar(std::string & read, std::pair<size_t, int> & mraPosLenReadF, std::string & ssr){
    this->ff = Sequence(read.substr(0, mraPosLenReadF.first));
    this->rf = Sequence(read.substr(mraPosLenReadF.first + mraPosLenReadF.second));
    this->mra = Sequence(read.substr(mraPosLenReadF.first, mraPosLenReadF.second));
    this->repuit = Sequence(ssr);
    this->effectiveSeq = Sequence(read);
    this->effectiveLen = read.length();
}

LocVar::LocVar(std::string & read, std::size_t & mraL, std::size_t & mraR, std::string & ssr){
    this->ff = Sequence(read.substr(0, mraL));
    this->rf = Sequence(read.substr(mraR));
    this->mra = Sequence(read.substr(mraL, mraR - mraL));
    this->effectiveSeq = Sequence(read);
    this->effectiveLen = read.length();
    this->repuit = Sequence(ssr);
}

LocVar::LocVar(std::string & name, std::string & ff, std::string & rf, std::string & repuit, std::string & mra){
    this->name = name;
    this->ff = Sequence(ff);
    this->rf = Sequence(rf);
    this->repuit = Sequence(repuit);
    this->mra = Sequence(mra);
    this->effectiveSeq = Sequence(ff + mra + rf);
    this->effectiveLen = ff.length() + mra.length() + rf.length();
}

void LocVar::print(){
    std::string msg = "name: " + name +  " ff: " + ff.mStr + " mra: " + mra.mStr + " rf: " + rf.mStr;
    cCout(msg, 'r');
}

void LocVar::printSnpMap(){
    std::string msg = "name: " + name +  " ff: " + ff.mStr + " mra: " + mra.mStr + " rf: " + rf.mStr;
    
    msg += "\nffSnpMap:\n";
    for(const auto & it : refSnpsSetffMap){
        for(const auto & it2 : it.second){
              msg += it.first + " : " + std::to_string(it2) + "\n";  
        }
    }
    
    msg += "rfSnpMap:\n";
    for(auto & it : refSnpsSetrfMap){
        for(const auto & it2 : it.second) {
                msg += it.first + " : " + std::to_string(it2) + "\n";
        }
    }
    cCout(msg, 'g');
}

LocVar::LocVar(std::string name, Sequence fp, Sequence rp, Sequence ff, Sequence rf, Sequence repuit, int nRep, unsigned int edCutoff, Sequence mra){
    this->name = name;
    this->fp = fp;
    this->rp = rp;
    this->ff = ff;
    this->rf = rf;
    this->repuit = repuit;
    this->nRep = nRep;
    this->edCutoff = edCutoff;
    this->mra = mra;
    this->effectiveSeq = ff.mStr + mra.mStr + rf.mStr;
    this->locSeq = fp.mStr + ff.mStr + mra.mStr + rf.mStr + rp.mStr;
    this->locLen = locSeq.length();
    this->effectiveLen = effectiveSeq.length();
    
}

Variance::Variance() {
    subMap.clear();
    delMap.clear();
    insMap.clear();
}

void Variance::cleanVar() {
    subMap.clear();
    delMap.clear();
    insMap.clear();
}


Genotype::Genotype() {
    numReads = 0;
    ssrF = "";
    ssrR = "";
}

void Genotype::print() {
    cCout("Num Reads: " + std::to_string(numReads), 'r');
    baseLocVar.print();
}

UnitedGenotype::UnitedGenotype(){
    marker = "";
    repuit = Sequence("");
    mra = Sequence("");
    effectiveSeq = Sequence("");
    //seqName = "";
    effectiveLen = 0;
    mraName = "";
    mraSVec.clear();
    genoSVec.clear();
    numReadsVec.clear();
    locVec.clear();
    effectiveSeqVec.clear();
    seqNameVec.clear();
}

LocSnp::LocSnp(){
    this->name = "";
    this->fp = Sequence("");
    this->rp = Sequence("");
    this->ref = Sequence("");
    this->snpPosSetHaplo.clear();
    this->snpPosSet.clear();
    this->refSnpPosSet.clear();
    this->numReads = 0;
    this->totReads = 0;
    this->totHaploReads = 0;
    this->snpsMap.clear();
    this->readsRatio = 0;
    this->genotype = "";
    this->puGeno = false;
    this->uGeno.preGenoMap.clear();
    this->uGeno.snpGenoMap.clear();
    this->haploVec.clear();
    this->isHaplotype = false;
    this->genoStr3 = "seqerr";
    this->ratioHaplo = 0;
}

void LocSnp::print(){
    std::string msg = "";
    msg = "name: " + name + " -> " + std::to_string(numReads) + "\n";
    msg += fp.mStr +  " :" + std::to_string(fp.mStr.length()) + "\n";
    msg += rp.mStr +  " :" + std::to_string(rp.mStr.length()) + "\n";
    msg += ref.mStr + " :" + std::to_string(ref.mStr.length()) + "\n";
    cCout(msg, 'r');
    msg = "";
    if (!snpsMap.empty()) {
        for (const auto & it : snpsMap) {
            msg += std::to_string(it.first) + " : " + it.second.first.mStr + "|" + it.second.second.mStr + "\n";
        }
    }
    cCout(msg, 'b');
}

std::string LocSnp::getGenotype(){
    if(snpsMap.empty()){
        return "ref";
    }
    std::string msg = "";
    for(const auto & it : snpsMap){
        msg += std::to_string(it.first) + "(" + it.second.first.mStr + "|" + it.second.second.mStr + ")";
    }
    return msg;
}

UnitedLocSnp::UnitedLocSnp(){
    this->preGenoMap.clear();
    this->snpGenoMap.clear();
    heter = false;
}

SeqVar::SeqVar(){
    this->seq = "";
    this->numReads = 0;
    this->indel = false;
    this->snpSet.clear();
}

LocSnp2::LocSnp2(){
    this->name = "";
    this->fp = Sequence("");
    this->rp = Sequence("");
    this->ref = Sequence("");
    this->snpPosSetHaplo.clear();
    this->snpPosSet.clear();
    this->refSnpPosSet.clear();
    //this->totPosSet.clear();
    this->totReads = 0;
    this->maxReads = 0;
    this->totHaploReads = 0;
    //this->totEffectReads = 0;
    this->ratioHaplo = 0;
    this->genoStr3 = "seqerr";
    this->isIndel = false;
    //this->haploVec.clear();
    //this->genoMap.clear();
    this->ft = Sequence("");
    this->rt = Sequence("");
    this->ssnpsMap.clear();
    this->baseErrorMap.clear();
    this->aveErrorRate = 0.000;
    this->seqVarVec.clear();
    this->status = std::make_pair(std::make_pair(false, false), false);
}

std::string LocSnp2::getHaploStr(bool snp2){
    std::string snpStr = "";
    if(ssnpsMap.empty()) return snpStr;
    if(snp2) {
        if (!status.first.second) {
            if (genoStr3 == "homo"){
                for (const auto &it : ssnpsMap)
                {
                    snpStr += it.second.snp1;
                }
            } else {
                for (const auto &it : ssnpsMap)
                {
                    snpStr += it.second.snp2;
                }
            }
                
        }
    } else {
        if(!status.first.first){
            for(const auto & it : ssnpsMap){
                snpStr += it.second.snp1;
            }
        }
    }
    return snpStr;
}

int LocSnp2::getNumSnps(){
    int num = 0;
    if (genoStr3 != "homo") {
        for (const auto & it : ssnpsMap) {
            if(seqVarVec.at(0).seq[it.first] != seqVarVec.at(1).seq[it.first]){
                num++;
            }
        }
    }
    return num;
}

std::string LocSnp2::getHaploStr(int index){
    std::string snpStr = "";
    if(seqVarVec.empty() || seqVarVec.at(index).indel || index >= seqVarVec.size()) return snpStr;
    for(const auto & it : seqVarVec.at(index).snpSet){
        snpStr += seqVarVec.at(index).seq[it];
    }
    return snpStr;
}

std::string LocSnp2::getSnpStr(bool snp2){
    std::string snpStr = "";
    if(ref.mStr.empty()) return snpStr;
    if(ssnpsMap.empty()) return snpStr;
    if(snp2) {
        if (!status.first.second && genoStr3 != "homo") {
            for (const auto & it : ssnpsMap) {
                snpStr += std::to_string(it.first + trimPos.first);
                snpStr += "(";
                snpStr += ref.mStr[it.first];
                snpStr += "|";
                snpStr += it.second.snp2;
                snpStr += ")";
            }
        }
    } else {
        if(!status.first.first){
            for(const auto & it : ssnpsMap){
                snpStr += std::to_string(it.first + trimPos.first);
                snpStr += "(";
                snpStr += ref.mStr[it.first];
                snpStr += "|";
                snpStr += it.second.snp1;
                snpStr += ")";
            }
        }
    }
    return snpStr;
}

std::string LocSnp2::getSnpStr(int index) {
    std::string snpStr = "";
    if (seqVarVec.empty() || seqVarVec.at(index).indel || index >= seqVarVec.size()) return snpStr;
    for (const auto & it : seqVarVec.at(index).snpSet) {
        snpStr += std::to_string(it + trimPos.first);
        snpStr += "(";
        snpStr += ref.mStr[it];
        snpStr += "|";
        snpStr += seqVarVec.at(index).seq[it];
        snpStr += ")";
    }
    return snpStr;
}

int LocSnp2::getHaploReads(bool haplo2){
    int num = 0;
    if(totHaploReads == 0) return num;
    if (genoStr3 == "homo") {
        if(haplo2){
            if (seqVarVec.size() > 1){
                num = seqVarVec.at(1).numReads;
            }
        } else {
            num = totHaploReads;
        }
    } else {
        if(haplo2){
            num = seqVarVec.at(1).numReads;
        } else {
            num = seqVarVec.at(0).numReads;
        }
    }
    return num;
}

int LocSnp2::getVarReads(int index){
    int num = 0;
    if(index <= seqVarVec.size()){
        num = seqVarVec.at(index).numReads;
    }
    return num;
}

double LocSnp2::getHaploReadsRatio(bool haplo2){
    if(haplo2){
        return getPer(getVarReads(1), getVarReads(0) + getVarReads(1));
    } else {
        if(seqVarVec.size() > 1){
            return getPer(getVarReads(0), getVarReads(0) + getVarReads(1));
        } else {
            return 1.0000;
        }
    }
}

double LocSnp2::getHaploReadsPer(bool haplo2){
    if(totHaploReads == 0 || totReads == 0) return 0.00;
    return getPer(getHaploReads(haplo2), totReads);
}

double LocSnp2::getReadsVarPer(int index){
    if(totReads == 0 || seqVarVec.at(index).numReads == 0) return 0.00;
    return getPer(seqVarVec.at(index).numReads, totReads);
}

std::string LocSnp2::getRatioStr() { 
    std::string str = ""; 
    //for(const auto & itv : ratioVec){
       // str += std::to_string(itv.first) + "|" + std::to_string(itv.second) + ";";
    //}
    if(!ratioVec.empty()){
        str += std::to_string(ratioVec.at(0).first) + "|" + std::to_string(ratioVec.at(0).second);
    }
    return str;
}

void LocSnp2::getBestRatio(){
    if (seqVarVec.size() > 2){
        bool go = false;
        if (status.second) {  // top 2 vairants have indels; only comparing the length, could have a bug
            go = true;
        } else {
            if(seqVarVec.at(0).seq.length() != seqVarVec.at(1).seq.length()){
                go = true;
            } else {
                std::set<int> tmpPosSet;
                tmpPosSet.insert(seqVarVec.at(0).snpSet.begin(), seqVarVec.at(0).snpSet.end());
                tmpPosSet.insert(seqVarVec.at(1).snpSet.begin(), seqVarVec.at(1).snpSet.end());
                for(const auto & it : tmpPosSet){
                    if(seqVarVec.at(0).seq.at(it) != seqVarVec.at(1).seq.at(it)){
                        ratioVec.push_back(getBaseFreqPair(it));
                    }
                }
                std::sort(ratioVec.begin(), ratioVec.end(),
                        [](const std::pair<int, double> &L, const std::pair<int, double> &R) {
                            return L.second < R.second;
                        });
            }
        }

        if(go){
            int len1 = seqVarVec.at(0).seq.length();
            int len2 = seqVarVec.at(1).seq.length();
            int len1n = 0;
            int len2n = 0;
            if (len1 != len2) {
                for(const auto & it : seqVarVec){
                    if(it.seq.length() == len1){
                        len1n++;
                    } else if (it.seq.length() == len2) {
                        len2n++;
                    }
                }
            }
            ratioVec.emplace_back(std::make_pair(0, getPer(std::max(len1n, len2n), (len1n + len2n), false)));
        }
    }
}

std::pair<int, double> LocSnp2::getBaseFreqPair(int pos){
    char c1n = seqVarVec.at(0).seq.at(pos);
    char c2n = seqVarVec.at(1).seq.at(pos);
    int n1 = 0;
    int n2 = 0;
    for (const auto &it : seqVarVec) {
        if(it.seq.at(pos) == c1n) {
            n1++;
        } else if (it.seq.at(pos) == c2n) {
            n2++;
        } else {

        }
    }
    return (std::make_pair((pos + trimPos.first), getPer(std::max(n1, n2), (n1 + n2), false)));
}

void LocSnp2::print(){
    std::string msg = "";
    msg = "name: " + name + " -> " + std::to_string(totReads) + "\n";
    msg += fp.mStr +  " :" + std::to_string(fp.mStr.length()) + "\n";
    msg += ft.mStr +  " :" + std::to_string(ft.mStr.length()) + "\n";
    msg += ref.mStr + " :" + std::to_string(ref.mStr.length()) + "\n";
    msg += rt.mStr +  " :" + std::to_string(rt.mStr.length()) + "\n";
    msg += rp.mStr +  " :" + std::to_string(rp.mStr.length()) + "\n";
    cCout(msg, 'r');
}

Sex::Sex(){
    this->sexMarker = "";
    this->primerF = Sequence("");
    this->primerR = Sequence("");
    this->refX = Sequence("");
    this->refY = Sequence("");
    this->maxReadsX = 0;
    this->maxReadsY = 0;
    this->totReadsX = 0;
    this->totReadsY = 0;
    this->minReadsSexVariant = 5;
    this->minReadsSexAllele = 10;
    this->haploRatio = 0.0;
    this->haplotype = false;
    this->haploIndel = false;
    this->haploStr = "inconclusive";
    this->baseErrorMapX.clear();
    this->baseErrorMapY.clear();
    this->aveErrorRateX = 0.0000;
    this->aveErrorRateY = 0.0000;
    this->mismatchesPF = 2;
    this->mismatchesPR = 2;
    this->mismatchesRX = 2;
    this->mismatchesRY = 2;
    this->lengthEqual = false;
    this->YXRatio = 0;
    this->YXRationCuttoff = 0.0001;
    this->sexMF = "Inconclusive";
    this->seqVarVecX.clear();
    this->seqVarVecY.clear();
}

void Sex::print(){
    std::cout << "name: " << sexMarker << "; primerF: " << primerF.mStr << "; primerR: " << primerR.mStr << 
            ";\n refx: " << refX.mStr << "; refy: " << refY.mStr << 
            ";\n readsX: " << maxReadsX << "; readsY: " << maxReadsY << "\n"; 
}

SeqVar Sex::getHaploVar(char sex, int index){
    SeqVar tmpSeqVar;
    if(sex == 'y' && !seqVarVecY.empty() && index < seqVarVecY.size()){
        tmpSeqVar =  seqVarVecY.at(index);
    } else if(sex == 'x' && !seqVarVecX.empty() && index < seqVarVecX.size()){
        tmpSeqVar =  seqVarVecX.at(index);
    }
    return tmpSeqVar;
}

std::string Sex::getSnpStr(char sex, int index){
    std::string snpStr = "";
    SeqVar tmp = getHaploVar(sex, index);
    
    if(tmp.snpSet.empty()){
        snpStr = "NA";
    } else {
        for (const auto & it : tmp.snpSet) {
            snpStr += std::to_string(it);
            snpStr += "(";
            snpStr += (sex == 'x' ? refX.mStr[it] : refY.mStr[it]);
            snpStr += "|";
            snpStr += tmp.seq[it];
            snpStr += ");";
            //            if(&it != &(*seqVarVecY.at(index).snpSet.rbegin())){
            //                snpStr += ";";
            //            }
        }
    }
    return snpStr;
}

std::string Sex::getFullRefX(){
    return std::string(primerF.mStr + refX.mStr + primerR.mStr);
}

std::string Sex::getFullRefY(){
    return std::string(primerF.mStr + refY.mStr + primerR.mStr);
}