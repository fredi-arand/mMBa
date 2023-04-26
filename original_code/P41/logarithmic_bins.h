#pragma once
#include <cstddef>

template <typename T>
struct LogarithmicBins{

    LogarithmicBins(T minValue, T maxValue, size_t bins)
      : minValue(minValue), maxValue(maxValue), bins(bins){}


    size_t to_bin(T const & value) const {
      long bin = (log(value/minValue)*T(bins))/log(maxValue/minValue);
      if(bin<0){ return 0;}
      if(bin>=bins){ return bins-1;}
      return bin;
    }

    T to_value(double bin) const {
      return minValue*pow(maxValue/minValue,bin/double(bins));
    }


    T const minValue, maxValue;
    size_t const bins;
};


template <typename T>
struct Bins{

    Bins(T minValue, T maxValue, size_t bins)
      : minValue(minValue), maxValue(maxValue), bins(bins){}


    size_t to_bin(T const & value) const {
      long bin = ((value-minValue)*T(bins))/(maxValue-minValue);
      if(bin<0){ return 0;}
      if(bin>=bins){ return bins-1;}
      return bin;
    }

    T to_value(double bin) const {
      return minValue+((maxValue-minValue)*bin)/double(bins);
    }


    T const minValue, maxValue;
    size_t const bins;
};
