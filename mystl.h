#pragma once

//
// mystl.h
// My Little Functions: Templates Are Magic
//
// This file contains small functions which just make life a bit easier
//

#include <map>
#include <string>
#include <sstream>

// Get std::map element. If element does not exist, return some default value
template <typename K, typename V>
V get_map_wd(const std::map <K,V> &m, const K &key, const V &defval ) {
   typename std::map<K,V>::const_iterator it = m.find( key );
   if ( it == m.end() ) {
      return defval;
   } else {
      return it->second;
   }
}

// Function to replace inclusions of oldStr in str with newStr
// You know, like 'string::replace' in normal languages
std::string string_replace(const std::string& str, const std::string& oldStr, const std::string& newStr);
// ToLower function
std::string lower(const std::string& str);

// Function to calculate gretest common divisor
// Not like it belongs here or anything...
// I haven't found a better place for it
// I'm really sorry
template <typename T>
T GCD(T a, T b) {
    T t;
    while (b != 0) {
        t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Convert anything to string
template <typename T> 
inline std::string any2str(const T& in) {
    std::stringstream s;
    s << in;
    return s.str();
}


template <typename T> 
inline T str2any(const char* in) {
    std::stringstream s(in);
    T out = T();
    s >> out;
    return out;
}

template <typename T>
inline T str2any(const std::string& in) {
	std::stringstream s(in);
	T out = T();
	s >> out;
	return out;
}


