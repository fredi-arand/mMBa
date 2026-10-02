#pragma once
//------------------------------------------------------------------------------
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
inline size_t numThreads() {
  return std::max(1u, std::thread::hardware_concurrency());
}
//------------------------------------------------------------------------------
// Calls f(i) for all i in [0, n) using all cores. Indices are handed out in
// small chunks, so uneven work per index is balanced between threads.
// f must be safe to call concurrently for different i.
template <typename F> void parallelFor(size_t n, F f) {
  if (n == 0)
    return;
  const size_t chunk = std::max<size_t>(1, n / (16 * numThreads()));
  std::atomic<size_t> next{0};
  auto worker = [&] {
    for (size_t first; (first = next.fetch_add(chunk)) < n;)
      for (size_t i = first; i < std::min(first + chunk, n); ++i)
        f(i);
  };

  std::vector<std::thread> threads(std::min(numThreads(), n) - 1);
  for (auto &thread : threads)
    thread = std::thread(worker);
  worker();
  for (auto &thread : threads)
    thread.join();
}
//------------------------------------------------------------------------------
// Sorts parts of [first, last) in parallel, then merges them pairwise.
// Same result as std::sort if cmp defines a total order.
template <typename It, typename Compare>
void parallelSort(It first, It last, Compare cmp) {
  const size_t n = last - first;
  const size_t parts = numThreads();
  if (parts == 1 || n < 10000) {
    std::sort(first, last, cmp);
    return;
  }

  std::vector<It> bounds;
  for (size_t i = 0; i <= parts; ++i)
    bounds.push_back(first + n * i / parts);
  parallelFor(parts,
              [&](size_t i) { std::sort(bounds[i], bounds[i + 1], cmp); });

  while (bounds.size() > 2) {
    parallelFor((bounds.size() - 1) / 2, [&](size_t i) {
      std::inplace_merge(bounds[2 * i], bounds[2 * i + 1], bounds[2 * i + 2],
                         cmp);
    });
    std::vector<It> merged;
    for (size_t i = 0; i < bounds.size(); i += 2)
      merged.push_back(bounds[i]);
    if (bounds.size() % 2 == 0)
      merged.push_back(bounds.back());
    bounds = merged;
  }
}
//------------------------------------------------------------------------------
} // namespace fred
