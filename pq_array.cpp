#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "pq_array.hpp"
#include "murmurhash.hpp"

using namespace std;

pq_array::pq_array (uint elements_per_queue, uint array_bits, uint c_bits, uint seed)
{
	this->seed = seed;
	this->c_bits = c_bits;
	this->array_bits = array_bits;
	this->queue_len = elements_per_queue;

	for (uint i = 0; i < (1U << array_bits); ++i)
	{
		this->array.push_back (std::vector<pq_element> ());
	}
}

void pq_array::add (uint64_t id, int32_t count)
{
	// Hash
	uint32_t hash = murmurhash (&id, this->seed);
	uint32_t array_idx = hash & ((1 << this->array_bits) - 1);
	//uint32_t int_id = (hash >> this->array_bits) & ((1UL << (72-c_bits))-1);
	uint32_t int_id = id;

	// Overwrite
	bool overwrite = false;
	for (size_t i = 0; i < this->array[array_idx].size (); ++i)
	{
		if (this->array[array_idx][i].id == int_id)
		{
			this->array[array_idx][i].count = count;
			overwrite = true;
		}
	}

	if (!overwrite)
		this->array[array_idx].push_back (pq_element (id, int_id, count));

	sort (this->array[array_idx].begin(), this->array[array_idx].end(),
			std::greater<pq_element>());

	if (this->array[array_idx].size () > queue_len)
		this->array[array_idx].pop_back ();
}

std::vector <uint32_t> pq_array::get_data ()
{
	std::vector <uint32_t> data;

	for (auto & v: array)
	{
		for (auto & e: v) data.push_back (e.count & ((1<<c_bits)-1) );
	}

	return data;
}

std::unordered_set <uint32_t> pq_array::get_id ()
{
	std::unordered_set <uint32_t> data;

	for (auto & v: array)
	{
		for (auto & e: v) data.insert (e.c_id);
	}

	return data;
}
