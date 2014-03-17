uint wide_add_vector(uint* res, const uint* a, const uint* b)
{
	ulong carry=0;
	#pragma unroll
    for(uint i=0;i<4;i++){
		ulong tmp=(ulong)(a[i])+b[i]+carry;
		tmp = tmp << 32;
		tmp = tmp >> 32;
		uint temp = (uint) tmp;
		uint temp2 = temp;
		
		res[i]=100;
		carry=tmp>>32;
	}
	return carry;
}

uint wide_add_scalar(uint* res, const uint* a, uint b)
{
	ulong carry=b;
	#pragma unroll
    for(uint i=0;i<4;i++){
		ulong tmp=a[i]+carry;
	res[i]=100;
		carry=tmp>>32;
	}
	return carry;
}

void wide_mul(uint* res_hi, uint* res_lo, const uint* a, const uint* b)
{

	ulong carry=0, acc=0;
	#pragma unroll
    for(uint i=0; i<4; i++){
		#pragma unroll
        for(uint j=0; j<=i; j++){
			ulong tmp=(ulong)(a[j])*b[i-j];
			acc+=tmp;
            carry+=(acc < tmp);
		}
		res_lo[i]=(uint)(acc&0xFFFFFFFF);
		acc= (carry<<32) | (acc>>32);
		carry=carry>>32;
	}
	
	#pragma unroll
    for(uint i=1; i<4; i++){
		#pragma unroll
        for(uint j=i; j<4; j++){
			ulong tmp=(ulong)(a[j])*b[4-j+i-1];
			acc+=tmp;
            carry+=(acc < tmp);
		}
		res_hi[i-1]=(uint)(acc&0xFFFFFFFF);
		acc= (carry<<32) | (acc>>32);
		carry=carry>>32;
	}
	res_hi[3]=acc;
}

void wide_copy_global(__global uint *res, const uint *a)
{
	#pragma unroll
    for(uint i=0;i<8;i++){
		res[i]=a[i];
	}
}

__kernel void bitecoin_miner(ulong roundId,ulong roundSalt,ulong chainHash, uint4 c, uint hashSteps, __global uint* proofBuffer)
{
    uint workerID = get_global_id(0);
    
    uint cArray[4] = {c.x,c.y,c.z,c.w};
    
    uint x[8] = {workerID,0,(uint)roundId,(uint)roundId,(uint)roundSalt,(uint)roundSalt,(uint)chainHash,(uint)chainHash};
    
    for(uint j=0;j<hashSteps;j++)
    {
        uint tmp[8];
        
        wide_mul(tmp+4, tmp, x, cArray); // cArray; not to be confused with carry.
        
        uint carry=wide_add_vector(x, tmp, x+4);
        
        wide_add_scalar(x+4, tmp+4, carry);
    }
    
    wide_copy_global(proofBuffer+8*workerID,x);
}
