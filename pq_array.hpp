#include <iostream>
#include <vector>
#include <unordered_set>

#ifndef PQ_ARRAY
#define PQ_ARRAY

class pq_element {
public:
	uint32_t c_id;
	uint32_t id;
	int32_t count;
	pq_element (uint32_t c_id, uint32_t id, int32_t count) : c_id (c_id), id (id), count (count) {};
	friend bool operator> (const pq_element& o1, const pq_element& o2)
	{
		return o1.count > o2.count;
	}
};

class pq_array
{
private:

	uint32_t seed;
	uint array_bits;
	uint c_bits;
	uint queue_len;
	std::vector<std::vector<pq_element>> array;
public:
	pq_array (uint elements_per_queue, uint array_bits, uint c_bits, uint seed);
	void add (uint64_t id, int32_t count);
	std::vector <uint32_t> get_data ();
	std::unordered_set <uint32_t> get_id ();
};

#endif // PQ_ARRAY
