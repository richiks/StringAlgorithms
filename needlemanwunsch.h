/**
  * @headerfile needlemanwunsch.h
  * @author Richik Vivek Sen (rsen9@gatech.edu)
  * @date 10/28/2019
  * @brief Implementation of Needleman Wunsch algorithm for string matching
  */

#ifndef NEEDLEMANWUNSCH_H
#define NEEDLEMANWUNSCH_H

#include <vector>
#include <algorithm> // For min
#include <iterator>  // For distance

/**
 * Function: NeedlemanWunschDistance(ForwardIterator1 begin1,
 *                                   ForwardIterator1 end1,
 *                                   ForwardIterator2 begin2,
 *                                   ForwardIterator2 end2);
 * Usage: cout << NeedlemanWunschDistance(str1.begin(), str1.end(),
 *                                        str2.begin(), str2.end());
 * ---------------------------------------------------------------------------
 * Given two sequences defined by [begin1, end1) and [begin2, end2), returns
 * the Needleman-Wunsch distance between those two sequences.
 */
template <typename ForwardIterator1, typename ForwardIterator2>
size_t NeedlemanWunschDistance(ForwardIterator1 begin1,
                               ForwardIterator1 end1,
                               ForwardIterator2 begin2,
                               ForwardIterator2 end2) {
  /* Begin by computing the sizes of the two ranges.  If we find that the
   * first range is smaller than the second range, exchange the two and
   * return that cost.
   */
  const size_t oneSize = size_t(std::distance(begin1, end1));
  const size_t twoSize = size_t(std::distance(begin2, end2));
  if (oneSize < twoSize)
    return NeedlemanWunschDistance(begin2, end2, begin1, end1);

  /* Construct a vector to hold the DP matching values, as well as a scratch
   * vector for use during each round.
   */
  std::vector<size_t> match(twoSize + 1), roundMatch(twoSize + 1);

  /* Base case: Cost of matching zero characters of the first string with some
   * number of characters of the second string is the number of characters in
   * the second string.
   */
  for (size_t i = 0; i < match.size(); ++i)
    match[i] = i;

  /* Inductive case: The cost of matching the first i characters of the first
   * string with the second string is defined above.
   */
  size_t i = 1;
  for (ForwardIterator1 itr1 = begin1; itr1 != end1; ++itr1, ++i) {
    /* The cost of matching the first i characters of the first string with
     * zero characters from the second string is i, since everything has to be
     * matched with a gap.
     */
    roundMatch[0] = i;

    /* Compute the recurrence. */
    size_t j = 1;
    for (ForwardIterator2 itr2 = begin2; itr2 != end2; ++itr2, ++j) {
      /* Compute the best we can do without applying a match. */
      size_t bestScore = 1 + std::min(roundMatch[j - 1], match[j]);

      /* If the characters match, update this to consider what happens when
       * we match them.
       */
      if (*itr1 == *itr2)
        bestScore = std::min(bestScore, match[j - 1]);

      /* Write the score out to the round vector. */
      roundMatch[j] = bestScore;
    }

    /* Update the resulting match score by swapping the scores from last round
     * and this round.  On the next round, we'll use the old match vector for
     * scratch space.
     */
    match.swap(roundMatch);
  }

  /* The final score is contained in the last slot of the round vector, which
   * corresponds to matching all characters of both strings.
   */
  return match.back();
}

#endif // NEEDLEMANWUNSCH_H
