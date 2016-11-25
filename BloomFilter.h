#pragma once

#include <stdint.h>
#include <vector>


struct BloomFilter {
public:
	BloomFilter(uint64_t size, uint8_t numHashes);

	void add(const char *data, std::size_t len);
	bool possiblyContains(const char *data, std::size_t len) const;

private:
	uint8_t m_numHashes;
	std::vector<bool> m_bits;
};
