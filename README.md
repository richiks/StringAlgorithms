# Rabin-Karp algorithm

 * An implementation of the Rabin-Karp string-matching algorithm.  Given two
 * strings - a text string T and a pattern string P, the Rabin-Karp algorithm
 * outputs the position of the first spot in which P matches T, or reports that
 * no such matching exists.
 *
 * The Rabin-Karp algorithm can be thought of as an improvement over the naive
 * string matching algorithm that tries to avoid unnecessary comparisons by
 * using a hash function.  In particular, suppose that we have an efficient way
 * (say, O(1)) for computing a hash function for each block of |P| consecutive
 * characters in the pattern string.  We can then get a very efficient string-
 * matching algorithm that works as follows.  For each block of |P| characters
 * in the string, compute its hash function in O(1).  If the hash value matches
 * the hash value of the pattern P, then check whether that particular sequence
 * matches the pattern P.  Assuming that we don't have too many hash 
 * collisions, this approach avoids a lot of unnecessary comparisons in which
 * a potential match is compared when there's no possible way for the substring
 * to match.
 *
 * The problem with this setup is that it makes a somewhat unrealistic
 * assumption - how can we compute hash functions for blocks of size O(|P|) in
 * time O(1)?  In general, this can't be done.  However, for specific choices
 * of hash functions, we can compute a hash function for a collection of blocks
 * of size O(|P|) in amortized O(1) time.  In particular, suppose that we have
 * a hash function with the following property: given a text string A0 A1 ...
 * An An+1 and given the hash code for A0 A1 ... An, then it's possible to
 * compute the hash code for the text string A1 A2 ... An An+1 in time O(1).
 * This type of hash function, called a "rolling hash function," gives us the
 * amortized guarantee we need.  We begin by computing the hash code for both
 * the pattern P and the first |P| characters of |T|, then scan over the text
 * of T, updating the rolling hash one character at a time and actually
 * checking for a match whenever the hash code of the current text string
 * matches the hash code of the pattern.
 *
 * The question now is - what's a good rolling hash funciton to use?  For this,
 * we'll use a hash function developed by Rabin and Karp specifically for this
 * purpose.  Suppose that the strings we're searching over are drawn from an
 * alphabet S = {0, 1, 2, ..., n - 1}.  Then we can treat any string in S* as
 * a number in base n.  For example, the string "12345" can be thought of as
 * number 12345, while the string "4144" would be the number 4144.  More
 * precisely, we can think of the string A1 A2 ... Ak as the base-n number
 *
 *                               k
 *                H[1 .. k] =   sum    Ai n^{k - i}
 *                             i = 1
 *
 * This then gives us a great way to update the hash code.  Given a text string
 * A1 ... Ak and a new character Ak+1, the hash code for the next value is
 *
 *                             k + 1
 *             H[2 .. k + 1] =  sum    Ai n^{k + 1 - i}
 *                             i = 2
 *
 * If we multiply the hash code for the previous string by n, we get
 *
 *                               k                     
 *               nH[1 .. k] =   sum    Ai n^{k + 1 - i}
 *                             i = 1                   
 *
 * And the difference of these two values is
 *
 *               H[1 .. k+1] - nH[0 .. k] = Ak+1 - n^k A1
 *
 * so H[2 .. k+1] = nH[1 .. k] + Ak+1 - n^k A1.  If we cache the value of
 * n^k, this can be computed in constant time.
 *
 * This is all great, but there's a catch - the hash code we're creating uses
 * O(lg k) bits to encode.  This isn't itself a problem, but on most modern
 * machines it can be a real hassle to work with variable-length hash codes. To
 * fix this, we modify the hash code so that all of the work is done modulo
 * some large prime q.  That is, the hash code is
 *
 *                               k
 *                H[1 .. k] =   sum    Ai n^{k - i}     (mod q)
 *                             i = 1
 *
 *                                k
 *                          = (  sum    (Ai n^{k - i}) mod q)  (mod q)
 *                              i = 1
 *
 * If we then work out the math once more about how to update the rolling hash
 * code, we get that
 *
 *                H[2 .. k+1] = nH[1 .. k] + Ak+1 - n^k A0 (mod q)
 *
 * Now, suppose that we pick a q such that nq fits into a single machine word
 * and such that n < q.  Then we know that nH[1 .. k] fits into a machine word,
 * since H[1 .. k] is modulo q.  We also know that Ak+1 fits into a machine
 * word, since it's a single character.  Finally, let's look at -n^k A1.  This
 * is equal to (((-n^k) mod q) A1) mod q, and also fits in a machine word
 * because A1 < n (since it's a character and is drawn from the set {0, 1, ..,
 * n-1}) and -n^k mod q < q.  Since all of these values fit into machine 
 * words, we can compute their sum in constant time.  In short, modding by q 
 * doesn't change the O(1) runtime for this step.  Interestingly, though, it 
 * does make this step much easier to implement.  Since nq fits into a word,
 * the hash code now takes up O(1) space, and more importantly space that can
 * be encoded using a single word rather than an array of words.
 *
 * The runtime of the Rabin-Karp algorithm is, in the worst case, O(|P||T|).
 * This happens if by sheer dumb luck we end up having a hash collision on
 * every length-|P| substring of T.  In the average case, though, it can be
 * shown that the number of hash collisions is expected O(1), in which case the
 * runtime of the algorithm is O(|P| + |T|), asymptotically much faster than
 * the naive implementation.  In practice, Rabin-Karp performs quite well on
 * most strings.

# Needleman-Wunsch algorithm

 * An implementation of the Needleman-Wunsch algorithm for optimal string
 * alignment.  The algorithm takes as input two strings, A and B, then 
 * computes the cost of an optimal alignment of the two strings formed by
 * inserting gaps at various points in the strings.  Gaps can be inserted
 * anywhere in the two strings, but all non-gap characters must align
 * perfectly with one another.  For example, to align TREE and THREE, we could
 * use any of the following alignments:
 *
 *                      T-REE   -T-REE-   -----TREE
 *                      THREE   T-HRE-E   THREE----
 *
 * Of these, the first is clearly the best, and the algorithm's goal is to
 * find and return it.
 *
 * The Needleman-Wunsch algorithm is interesting for a variety of reasons.
 * Historically, it was one of the first major biological algorithms to use
 * dynamic programming, and is often one of the first DP algortihms taught in
 * most algorithms classes.  It is also interesting because a naive
 * implementation takes O(mn) memory, while a more clever version can use
 * O(min{m, n}) memory (this version is implemented here).
 *
 * The idea behind the algorithm is, like the Levenshtein distance algorithm,
 * to consider the optimal way of matching the first characters of the
 * strings.  There are three possible options:
 *
 * 1. If the first two characters match, the cost is the cost of matching
 *    the rest of the string.
 * 2. Otherwise, we could match the first character of the first string with
 *    a gap, then match the remainder of that string with the second string.
 *    The cost is one plus the cost of that matching.
 * 3. Finally, we could match the first character of the second string with
 *    a gap, for a total cost of one plus the cost of matching the rest of the
 *    second string with the first string.
 *
 * If we let C(i, j) be the cost of matching the first m characters of the
 * first string with the first n characters of the second string, this
 * recurrence relation is formalized as
 *
 * C(0, j) = j     (The only way to match one string with an empty string is
 *                  to match the rest of the string with blanks)
 * C(i, 0) = i     (Same)
 * C(i, j) = min{ C(i - 1, j - 1)  (if A[i] == B[j], one-indexed),
 *                1 + C(i - 1, j),
 *                1 + C(i, j - 1) }
 *
 * Because this recurrence displays a large degree of overlapping subproblems,
 * it's a perfect candidate for a dynamic programming solution.  We initialize
 * a grid of size (m + 1) x (n + 1) to zero, initialize the base cases (for 
 * when i == 0), and then proceed upward and across filling it in.  This uses 
 * O(mn) memory and runs in O(mn) time.
 *
 * However, we can do much better than this.  Look at the effect of the
 * recurrence as a function of n.  Each term touched by the recurrence looks
 * only at lower values of m in the same row as n, or at the previous row of
 * n.  This means that we don't actually need to store the entire table, and
 * in fact only need access to the last row.  This uses memory O(n), which is
 * significantly better than before.
 *
 * But we can do even better than this!  Let m and n be the lengths of strings
 * A and B, respectively.  If m < n, then we can exchange strings A and B and
 * then use only O(m) memory.  This means that our memory usage is thus
 * O(min{m, n}), much better than what we started with.
 *
 * Some variants of this problem assign different weights to matchings of
 * various characters with blanks, or allows each character to match against
 * each other character for some penalty.  We don't consider this here, but
 * it's quite easy to adapt the algorithm to handle this.
