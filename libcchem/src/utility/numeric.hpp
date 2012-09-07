#ifndef CCHEM_UTILITY_NUMERIC_HPP
#define CCHEM_UTILITY_NUMERIC_HPP

namespace utility {

  template<class T> static inline T binomial2(T t) {
    return (t * t + t) / 2;
  }


  template<class T> static inline T pow2(T t) {
    return t*t;
  }

  template<class T> static inline T pow4(T t) {
    return pow2(t*t);
  }

  template<typename T> static inline T factorial2(T t) {
    if (t < 2) return t;
    return t * factorial2(t - 2);
  }

  template<typename T> static inline void permute(int i, int j, int k, std::vector<T> &v) {
    if (i == j && j == k) {
      v.push_back(T(i, j, k));
    }
    else if (j == k) {
      v.push_back(T(i, j, k));
      v.push_back(T(j, i, k));
      v.push_back(T(k, j, i));
    }
    else if (i == j) {
      v.push_back(T(i, j, k));
      v.push_back(T(i, k, j));
      v.push_back(T(k, j, i));
    }
    else {
      v.push_back(T(i, j, k));
      v.push_back(T(i, k, j));
      v.push_back(T(j, i, k));
      v.push_back(T(k, i, j));
      v.push_back(T(j, k, i));
      v.push_back(T(k, j, i));
    }
  }

}

#endif
