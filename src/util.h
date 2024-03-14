#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>
#include <sstream>
#include <array>
#include <ctype.h>
#include <stdio.h>
#include <map>
#include <cmath>
#include "common.h"

using namespace std;

extern mutex logmtx;

inline void loginfo(const string s) {
    logmtx.lock();
    time_t tt = time(NULL);
    tm* t = localtime(&tt);
    //cerr<<"["<<t->tm_hour<<":"<<t->tm_min<<":"<<t->tm_sec<<"] "<<s<<endl;
    //fprintf(stderr, "[\033[1;35m%02d:%02d:%02d\033[0m] %s\n", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    fprintf(stderr, "[%02d:%02d:%02d] %s\n", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    logmtx.unlock();
}

template <typename T>
void cCout(const T & str, char color = 'd', bool newLine = true) {
    logmtx.lock();
    if (color == 'r') {
        if (newLine) {
            std::cout << "\033[1;31m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;31m" << str << "\033[0m";
        }
    } else if (color == 'g') {
        if (newLine) {
            std::cout << "\033[1;32m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;32m" << str << "\033[0m";
        }

    } else if (color == 'y') {
        if (newLine) {
            std::cout << "\033[1;33m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;33m" << str << "\033[0m";
        }

    } else if (color == 'b') {
        if (newLine) {
            std::cout << "\033[1;34m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;34m" << str << "\033[0m";
        }

    } else if (color == 'm') {
        if (newLine) {
            std::cout << "\033[1;35m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;35m" << str << "\033[0m";
        }

    } else if (color == 'c') {
        if (newLine) {
            std::cout << "\033[1;36m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;36m" << str << "\033[0m";
        }

    } else if (color == 'w') {
        if (newLine) {
            std::cout << "\033[1;37m" << str << "\033[0m" << std::endl;
        } else {
            std::cout << "\033[1;37m" << str << "\033[0m";
        }

    } else if (color = 'd') {
        if (newLine) {
            std::cout << str << std::endl;
        } else {
            std::cout << str;
        }
        //        if (newLine) {
        //            std::cout << "\033[1;30m" << str << "\033[0m" << std::endl;
        //        } else {
        //            std::cout << "\033[1;30m" << str << "\033[0m";
        //        }
    } else {

    }
    logmtx.unlock();
}

template <typename T1, typename T2>
void cCout(const T1 & str1, const T2 & str2, char color = 'r') {
    logmtx.lock();
    if (color == 'r') {
        std::cout << "\033[1;31m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'g') {
        std::cout << "\033[1;32m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'y') {
        std::cout << "\033[1;33m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'b') {
        std::cout << "\033[1;34m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'm') {
        std::cout << "\033[1;35m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'c') {
        std::cout << "\033[1;36m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'w') {
        std::cout << "\033[1;37m" << str1 << " -> " << str2 << "\033[0m\n";
    } else {
        std::cout << "\033[1;30m" << str1 << " -> " << str2 << "\033[0m\n";
    }
    logmtx.unlock();
}

template <typename T>
std::string addColor(const T & str, char color = 'd', bool newLine = true) {
    std::stringstream ss;
    if (color == 'r') {
        if (newLine) {
            ss << "\033[1;31m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;31m" << str << "\033[0m";
        }
    } else if (color == 'g') {
        if (newLine) {
            ss << "\033[1;32m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;32m" << str << "\033[0m";
        }

    } else if (color == 'y') {
        if (newLine) {
            ss << "\033[1;33m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;33m" << str << "\033[0m";
        }

    } else if (color == 'b') {
        if (newLine) {
            ss << "\033[1;34m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;34m" << str << "\033[0m";
        }

    } else if (color == 'm') {
        if (newLine) {
            ss << "\033[1;35m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;35m" << str << "\033[0m";
        }

    } else if (color == 'c') {
        if (newLine) {
            ss << "\033[1;36m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;36m" << str << "\033[0m";
        }

    } else if (color == 'w') {
        if (newLine) {
            ss << "\033[1;37m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;37m" << str << "\033[0m";
        }

    } else if (color = 'd') {
        if (newLine) {
            ss << "\033[1;30m" << str << "\033[0m\n";
        } else {
            ss << "\033[1;30m" << str << "\033[0m";
        }
    } else {

    }

    return ss.str();
}

inline void error_exit(const string& msg) {
    cerr << "\033[1;31mERROR: " << msg << "\033[0m" << endl;
    exit(-1);
}

inline char complement(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with(string const & value, string const & starting) {
    if (starting.size() > value.size()) return false;
    return equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with(string const & value, string const & ending) {
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string& str) {
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos) {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline string getPrefix(std::string & str) {
    for (const auto & it : fileTypes) {
        if (ends_with(str, it)) {
            return str.substr(0, str.length() - it.length());
        }
    }
    return str;
}

inline int split(const string& str, vector<string>& ret_, string sep = ",") {
    if (str.empty()) {
        return 0;
    }

    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;

    while (pos_begin != string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        ret_.push_back(tmp);
        tmp.clear();
    }
    return 0;
}

inline string replace(const string& str, const string& src, const string& dest) {
    string ret;

    string::size_type pos_begin = 0;
    string::size_type pos = str.find(src);
    while (pos != string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string reverse(const string& str) {
    string ret(str.length(), 0);
    for (int pos = 0; pos != str.length(); ++pos) {
        ret[pos] = str[str.length() - pos - 1];
    }
    return ret;
}

inline string basename(const string& filename) {
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;
    else if (pos == filename.length() - 1)
        return ""; // a bad filename
    else
        return filename.substr(pos + 1, filename.length() - pos - 1);
}

inline string dirname(const string& filename) {
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos) {
        return "./";
    } else
        return filename.substr(0, pos + 1);
}

inline string joinpath(const string& dirname, const string& basename) {
    if (dirname[dirname.length() - 1] == '/') {
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

//Check if a string is a file or directory

inline bool file_exists(const string& s) {
    bool exists = false;
    if (s.length() > 0) {
        struct stat status;
        int result = stat(s.c_str(), &status);
        if (result == 0) {
            exists = true;
        }
    }
    return exists;
}


// check if a string is a directory

inline bool is_directory(const string& path) {
    bool isdir = false;
    struct stat status;
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat(path.c_str(), &status);
    if (status.st_mode & S_IFDIR) {
        isdir = true;
    }
    // #endif
    return isdir;
}

inline void check_file_valid(const string& s) {
    if (!file_exists(s)) {
        cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
        exit(-1);
    }
    if (is_directory(s)) {
        cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
        exit(-1);
    }
}

inline void check_file_writable(const string& s) {
    string dir = dirname(s);
    if (!file_exists(dir)) {
        cerr << "ERROR: '" << dir << " doesn't exist. Create this folder and run this command again." << endl;
        exit(-1);
    }
    if (is_directory(s)) {
        cerr << "ERROR: '" << s << "' is not a writable file, quit now" << endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string

inline string str_keep_alpha(const string& s) {
    string new_str;
    for (size_t it = 0; it < s.size(); it++) {
        if (isalpha(s[it])) {
            new_str += s[it];
        }
    }
    return new_str;
}


// Remove invalid sequence characters from a string

inline void str_keep_valid_sequence(string& s, bool forceUpperCase = false) {
    size_t total = 0;
    const char case_gap = 'a' - 'A';
    for (size_t it = 0; it < s.size(); it++) {
        char c = s[it];
        if (forceUpperCase && c >= 'a' && c <= 'z') {
            c -= case_gap;
        }
        if (isalpha(c) || c == '-' || c == '*') {
            s[total] = c;
            total++;
        }
    }

    s.resize(total);
}

inline int find_with_right_pos(const string& str, const string& pattern, int start = 0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline void str2upper(string& s) {
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper);
}

inline void str2lower(string& s) {
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower);
}

inline char num2qual(int num) {
    if (num > 127 - 33)
        num = 127 - 33;
    if (num < 0)
        num = 0;

    char c = num + 33;
    return c;
}

inline int extractIntegerWords(string & str) {
    stringstream ss;
    ss << str;

    /* Running loop till the end of the stream */
    string temp;
    int max = 0;
    int found;
    while (!ss.eof()) {

        /* extracting word by word from stream */
        ss >> temp;

        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found) {
            //cout << found << "\n";
            if (found > max) {
                max = found;
            }
        }
        /* To save from space at the end of string */
        temp = "";
    }

    return max;
}

inline void splitStr(const string & str, std::vector<std::string> & ret_, string sep = "\t") {
    if (str.empty()) return;

    std::string tmp;
    std::string::size_type pos_begin = str.find_first_not_of(sep);
    std::string::size_type sep_pos = 0;

    while (pos_begin != std::string::npos) {
        sep_pos = str.find(sep, pos_begin);
        if (sep_pos != std::string::npos) {
            tmp = str.substr(pos_begin, sep_pos - pos_begin);
            pos_begin = sep_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = sep_pos;
        }
        ret_.emplace_back(tmp);
        tmp.clear();
    }
    ret_.shrink_to_fit();

    return;
}

inline int findDigitLen(std::string str) {
    int len = 0;
    for (int i = 0; i != str.length(); ++i) {
        if (isdigit(str[i])) {
            len++;
        } else {
            break;
        }
    }
    return len;
}

inline std::string findSingSSR(std::string & cssr) {
    std::string newStr;
    int pos = 0;
    while (pos != cssr.length()) {
        if (cssr[pos] == '(') {
            auto len = findDigitLen(cssr.substr(pos + 5));
            newStr += cssr.substr(pos, len + 5);
            pos += (5 + len);
        } else {

            if (pos != cssr.length() - 1) {
                if (cssr[pos] == cssr[pos + 1]) {
                    newStr += "(";
                    newStr += cssr[pos];
                    int rep = 1;
                    while (pos != cssr.length() - 1) {
                        if (cssr[pos] == cssr[pos + 1]) {
                            rep++;
                            pos++;
                            if (pos == cssr.length() - 1) {
                                newStr += ")";
                                newStr += std::to_string(rep);
                                return (newStr);
                            }
                        } else {
                            pos++;
                            newStr += ")";
                            newStr += std::to_string(rep);
                            break;
                        }
                    }

                } else {
                    newStr += cssr[pos];
                    pos++;
                }

            } else {
                newStr += cssr[pos];
                pos++;
            }
        }

    }

    return newStr;
}

inline std::string getGenotype(std::string & mra, std::string & ssr, bool adp = false) {
    std::string genStr = "";
    size_t startPos = mra.find(ssr);
    size_t endPos = startPos;
    int count = 0;
    if (startPos != 0) {
        genStr.append(mra.substr(0, startPos));
    }

    while (startPos != std::string::npos) {
        count++;
        endPos = mra.find(ssr, startPos + ssr.length());
        if (endPos != std::string::npos) {
            if (endPos - startPos == ssr.length()) {

            } else {
                if (adp) {
                    genStr.append("|(" + ssr + ")" + std::to_string(count) + "|" + mra.substr(startPos + ssr.length(), endPos - (startPos + ssr.length())));
                } else {
                    genStr.append("(" + ssr + ")" + std::to_string(count) + mra.substr(startPos + ssr.length(), endPos - (startPos + ssr.length())));
                }
                count = 0;
            }
        } else {
            if (adp) {
                genStr.append("|(" + ssr + ")" + std::to_string(count) + "|" + mra.substr(startPos + ssr.length()));
            } else {
                genStr.append("(" + ssr + ")" + std::to_string(count) + mra.substr(startPos + ssr.length()));
            }
            count = 0;
        }
        startPos = endPos;
    }

    //if(adp){
    return genStr;
    //} else {
    //return findSingSSR(genStr);
    //}

}

inline std::string getGenotype(std::string & mra, std::string & ssr, std::string & ssr2) {
    if (ssr2.empty()) {
        return getGenotype(mra, ssr);
    } else {
        std::string genStr = getGenotype(mra, ssr, true);
        std::vector<std::string> splitGenStrVec;
        splitStr(genStr, splitGenStrVec, "|");
        genStr = "";
        for (auto & it : splitGenStrVec) {
            if (isalpha(it[0])) {
                genStr.append(getGenotype(it, ssr2));
            } else {
                genStr.append(it);
            }
        }
        return (genStr.empty() ? mra : genStr);
    }
}

inline int getMaxNumRep(string & str) {
    int max = 0;

    std::string tmpNum = "";
    size_t startPos = str.find(")");
    size_t endPos;

    while (startPos != std::string::npos) {
        auto it = str[startPos + 1];
        if (isdigit(it)) {
            tmpNum.push_back(it);
            startPos++;
            if (startPos == str.size()) {
                break;
            } else {
                continue;
            }
        } else {
            if (max < stoi(tmpNum)) {
                max = stoi(tmpNum);
            }
            tmpNum = "";
            startPos = str.find(")", startPos);
        }
    }

    return max;
}

inline std::string getMraBase(std::string & genStr) {
    std::stringstream ss;
    for (int i = 0; i != genStr.length(); ++i) {
        auto tmp = genStr[i];
        if (!isdigit(tmp)) {
            ss << tmp;
        }
    }
    return ss.str();
}

inline string trimStr(const string& str) {
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos) {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

template<typename K, typename V>
std::pair<K, V> getMinKeyValue(const std::map<K, V> & map, bool reverse = false) {
    return *std::min_element(map.begin(), map.end(),
            [](std::pair<K, V> const & x, std::pair<K, V> const & y) {
                return x.second < y.second;
            });
}

template<typename K, typename V>
std::pair<K, V> getMaxKeyValue(const std::map<K, V> & map, bool reverse = false) {
    if (reverse) {
        return *std::max_element(map.rbegin(), map.rend(),
                [](std::pair<K, V> const & x, std::pair<K, V> const & y) {
                    return x.second < y.second;
                });
    } else {
        return *std::max_element(map.begin(), map.end(),
                [](std::pair<K, V> const & x, std::pair<K, V> const & y) {
                    return x.second <= y.second;
                });
    }
}

template<typename K, typename V>
std::vector<std::pair<K, V>> getTop2MaxKeyValueVec(const std::map<K, V>& map) {
    std::map<K, V> tmp = map;
    std::pair<K, V> m1 = getMaxKeyValue(tmp);
    tmp.erase(m1.first);
    std::vector<std::pair<K, V>> tmpVec;

    tmpVec.push_back(m1);

    if (!tmp.empty()) {
        std::pair<K, V> m2 = getMaxKeyValue(tmp, true);
        tmp.clear();
        tmpVec.push_back(m2);
    }

    return (tmpVec);
}

template<typename K, typename V>
std::map<K, V> getTop2MaxKeyValue(const std::map<K, V> & map) {
    std::map<K, V> tmp = map;
    std::pair<K, V> m1 = getMaxKeyValue(tmp);
    tmp.erase(m1.first);

    if (tmp.empty()) {
        tmp.clear();
        tmp = {m1};
    } else {
        std::pair<K, V> m2 = getMaxKeyValue(tmp, true);
        tmp.clear();
        tmp = {m1, m2};
    }

    return (tmp);
}

template<typename K, typename V>
std::map<K, V> get2Peaks(const std::map<K, V> & map, double hlRatio, double hlRatio2) {
    auto it = map.begin();
    auto nReads = it->second;
    auto nReadsPeak = nReads;
    auto loc = it->first;
    auto locPeak = loc;
    std::map<K, V> res;
    int i = 1;
    while (i < map.size()) {
        i++;
        std::advance(it, 1);
        if (nReadsPeak < it->second) {
            nReadsPeak = it->second;
            locPeak = it->first;
        }
        if (it->second >= nReads) {
            nReads = it->second;
            loc = it->first;
        } else {
            res[loc] = nReads;
            nReads = it->second;
            loc = it->first;
            break;
        }
    }
    if (res.empty()) {
        res[locPeak] = nReadsPeak;
        if (nReadsPeak == nReads && locPeak != loc) {
            res[loc] = nReads;
        }
        return (res);
    }

    bool cont = false;
    if (i == map.size()) {
        if ((double) nReads / nReadsPeak >= hlRatio) {
            res[loc] = nReads;
            return (res);
        } else {
            cont = true;
        }

    } else if (i == map.size() - 1) {
        i++;
        std::advance(it, 1);
        if ((double) nReads / nReadsPeak >= hlRatio && (double) (it->second) / nReadsPeak >= hlRatio2) {
            if (nReads > it->second) {
                res[loc] = nReads;
            } else {
                res[it->first] = it->second;
            }

            return (res);
        } else {
            cont = true;
        }
    } else {
        i++;
        std::advance(it, 1);
        cont = true;
    }

    if (cont) {
        nReads = it->second;
        loc = it->first;
        while (i < map.size()) {
            i++;
            std::advance(it, 1);
            if (it->second >= nReads) {
                nReads = it->second;
                loc = it->first;
            }
        }

        std::map<K, V> map2;
        for (const auto & itt : map) {
            if (itt.first >= loc) {
                map2[itt.first] = itt.second;
            }
        }
        res.insert(getMaxKeyValue(map2, true));
    }
    return (res);
}

template<typename K, typename V>
double getPer(K k, V v, bool per = true){
    if(v == 0) return 0.00;
    double res = 0.0;
    if(per){
        res = std::round((static_cast<double>(k * 100.0) / v) * 10000.0) / 10000.0;
    } else {
        res = std::round((static_cast<double>(k) / v) * 10000.0) / 10000.0;
    }
    return res;
}
#endif /* UTIL_H */
