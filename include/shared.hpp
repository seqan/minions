#pragma once

//
/*! \brief Function that ensures random hashes, based on https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
 *  \param hash_value The hash_value that should be transformed.
 *  \param seed       The seed.
 */
uint64_t fnv_hash(uint64_t hash_value, uint64_t seed)
{
    // If seed is 0, then the hash value is just returned.
    if (seed == 0)
        return hash_value;

    constexpr static uint64_t default_offset_basis = 0xcbf29ce484222325;
    constexpr static uint64_t prime                = 0x100000001b3;

    uint64_t hashed = hash_value;
    std::ostringstream os;
    os << hash_value;
    std::string oss = os.str();

    for (int i = 0; i < oss.size(); i++)
    {
        hashed = hashed * prime;
        hashed= hashed ^ oss[i];
    }

    return hashed;
}
