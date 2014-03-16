uint min2(uint a,uint b)
{
    uint mask = 0xC0000000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    for(int i=0;i<16;i++)
    {
        x = a & mask;
        y = b & mask;
        
        // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
        
        output |= y ^ ((x ^ y) & -(x < y));
        
        mask >>= 2;
    }
    
    return output;
}

uint max2(uint a,uint b)
{
    uint mask = 0xC0000000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    for(int i=0;i<16;i++)
    {
        x = a & mask;
        y = b & mask;
        
        // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
        
        output |= x ^ ((x ^ y) & -(x < y));
        
        mask >>= 2;
    }
    return output;
}

uint min4(uint a,uint b)
{
    uint mask = 0xF0000000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    for(int i=0;i<8;i++)
    {
        x = a & mask;
        y = b & mask;
        
        // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
        
        output |= y ^ ((x ^ y) & -(x < y));
        
        mask >>= 4;
    }
    
    return output;
}

uint max4(uint a,uint b)
{
    uint mask = 0xF0000000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    for(int i=0;i<8;i++)
    {
        x = a & mask;
        y = b & mask;
        
        // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
        
        output |= x ^ ((x ^ y) & -(x < y));
        
        mask >>= 4;
    }
    
    return output;
}

uint min8(uint a,uint b)
{
    uint mask = 0xFF000000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    for(int i=0;i<4;i++)
    {
        x = a & mask;
        y = b & mask;
        
        // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
        
        output |= y ^ ((x ^ y) & -(x < y));
        
        mask >>= 8;
    }
    
    return output;
}

uint max8(uint a,uint b)
{
    uint mask = 0xFF000000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    for(int i=0;i<4;i++)
    {
        x = a & mask;
        y = b & mask;
        
        // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
        
        output |= x ^ ((x ^ y) & -(x < y));
        
        mask >>= 8;
    }
    
    return output;
}

uint min16(uint a,uint b)
{
    uint mask = 0xFFFF0000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    x = a & mask;
    y = b & mask;
    
    output |= y ^ ((x ^ y) & -(x < y));
    
    mask >>= 16;
    
    x = a & mask;
    y = b & mask;
    
    output |= y ^ ((x ^ y) & -(x < y));
    
    // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
    
    return output;
}

uint max16(uint a,uint b)
{
    uint mask = 0xFFFF0000;
    uint output = 0x00000000;
    
    uint x;
    uint y;
    
    x = a & mask;
    y = b & mask;
    
    output |= y ^ ((x ^ y) & -(x < y));
    
    mask >>= 16;
    
    x = a & mask;
    y = b & mask;
    
    output |= x ^ ((x ^ y) & -(x < y));
    
    // Non-branching minimum taken from http://graphics.stanford.edu/~seander/bithacks.html
    
    return output;
}

__kernel void bitecoin_miner(ulong roundId,ulong roundSalt,ulong chainHash,ulong cHi,ulong cLo,uint maxIndices, uint hashSteps,__global uint* indexBuffer, __global ulong* proofBuffer)
{
    uint workerID = get_global_id(0);
    
    // Build random index of indexsize
    // curr=curr+1+(rand()%10) for 1 to maxIndices; how do we do this dynamically??
    
    
    
}


