#pragma once

#include <vector>
#include <random>
#include <stddef.h>

template <typename T>
struct Histogram{

    Histogram(T minValue, T maxValue, size_t bins){
      this->minValue = minValue;
      this->maxValue = maxValue;
      this->bins = bins;
      counts.resize(bins,0);
    }

    size_t & operator()(T const & value){
      long bin = to_bin(value);
      if(bin<0){ return counts[0]; }
      if(bin>=bins){ return counts[bins-1]; }
      return counts[bin];
    }

    size_t operator()(T const & value) const {
      long bin = to_bin(value);
      if(bin<0){ return counts[0]; }
      if(bin>=bins){ return counts[bins-1]; }
      return counts[bin];
    }

    T to_value(int const & bin) const {
      return (float(bin*float(maxValue-minValue))/float(bins)+float(minValue));
    }

    long to_bin( T const & value ){
      return (T(bins)*(value-minValue))/(maxValue-minValue);
    }

    T minValue;
    T maxValue;
    size_t bins;

    std::vector<size_t> counts;
};

//*************************************************************************************************

Histogram<double> example_histogram()
{

  using namespace std;

  Histogram<double> histogram(-0.5,255.5,256);

  size_t samples = pow(256,3);
  default_random_engine generator;
  normal_distribution<double> distribution0(80.0,10.0);

  size_t sample = 0;
  while(sample != samples)
  {
    double value = distribution0(generator);
    if(value<-0.5 || value > 255.5)
      continue;

    ++sample;
    ++histogram(value);
  }

  normal_distribution<double> distribution1(120.0,10.0);
  sample = 0;
  while(sample!= samples/100)
  {
    double value = distribution1(generator);
    if(value<-0.5 || value > 255.5)
      continue;

    ++sample;
    ++histogram(value);
  }

  return histogram;
}
