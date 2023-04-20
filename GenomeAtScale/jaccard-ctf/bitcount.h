// Trying different implementations of popcount
// Can use hardware instruction based on the architecture
// Source: Stack Overflow, Wiki, etc

template <typename T>
uint32_t popcount(T x)
{
//throw
  throw std::runtime_error("bitmask type not recognised");
}

template <>
uint32_t popcount(uint32_t i)
{
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

template<>
uint32_t popcount(uint64_t i)
{
  i = i - ((i >> 1) & 0x5555555555555555);
  i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
  i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
  i = i + (i >> 8);
  i = i + (i >> 16);
  i = i + (i >> 32);
  return (uint32_t)i & 0x7f;
}

/*
 * lookup table implementation
 * bit faster, as tested on Stampede2 with m=1M, n=10k, 16 nodes (~19% faster)
unsigned char wordbits[65536] = {bitcounts of ints between 0 and 65535};
uint32_t popcount(uint32_t i)
{
  return(wordbits[i&0xFFFF] + wordbits[i>>16]);
}
*/
