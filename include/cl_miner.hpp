// Header files for windows compilation
#ifdef _WIN32
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#endif

// Shared Headers
#define __CL_ENABLE_EXCEPTIONS

#include <stdexcept>
#include <cmath>
#include <stdint.h>
#include <memory>
#include <cstdio>
#include <fstream>
#include <streambuf>
#include <iostream>
#include <tuple>
#include <sys/time.h>

#include "bitecoin_hashing.hpp"


//#include "tbb/parallel_for.h"
//#include "tbb/task_group.h"

#include "CL/cl.hpp"
#include "cl_utilities.hpp"

bitecoin::bigint_t do_gpu_mining(uint32_t microrounds, const bitecoin::Packet_ServerBeginRound *pParams, uint64_t chainHash, cl_instance instance)
{
    bitecoin::bigint_t dummy;
    wide_zero(8,dummy.limbs);
    return dummy;
}

//void process_opencl_packed_line(int levels, unsigned w, unsigned bits,std::vector<uint32_t>& gpuReadOffsets, std::vector<uint32_t>& gpuWriteOffsets, uint32_t* pixelsIn, uint32_t* pixelsOut,std::vector<uint32_t> aboveOverrides,std::vector<uint32_t> belowOverrides,std::tuple<cl::Kernel,cl::Kernel,std::vector<cl::Buffer*>,cl::CommandQueue,cl::NDRange,cl::NDRange,cl::NDRange> cl_instance)
//{
//    size_t cbBuffer=(w*bits)/8;
//    
//    cl::Kernel erodeKernel = std::get<0>(cl_instance);
//    cl::Kernel dilateKernel = std::get<1>(cl_instance);
//    std::vector<cl::Buffer*> gpuBuffers = std::get<2>(cl_instance);
//    cl::CommandQueue queue = std::get<3>(cl_instance);
//    cl::NDRange offset = std::get<4>(cl_instance);
//    cl::NDRange globalSize = std::get<5>(cl_instance);
//    cl::NDRange localSize = std::get<6>(cl_instance);
//    
//    queue.enqueueWriteBuffer(*gpuBuffers[0], CL_FALSE, cbBuffer*gpuWriteOffsets[0], cbBuffer, pixelsIn);
//    
//    dilateKernel.setArg(5, *gpuBuffers[i+1]);
//    
//    queue.enqueueNDRangeKernel(dilateKernel, offset, globalSize, localSize);
//    
//    queue.enqueueReadBuffer(*gpuBuffers[2*abs(levels)], CL_TRUE, cbBuffer*gpuReadOffsets[2*abs(levels)], cbBuffer, pixelsOut);
//
//    queue.enqueueBarrier();
//}

//void transform(int deviceNumber, int levels, unsigned w, unsigned h, unsigned bits)
//{
//    uint64_t cbinput=uint64_t(w)*bits/8;
//    
//    
//    std::vector<uint64_t> input(2*(cbinput/8));
//    std::vector<uint64_t> output(2*(cbinput/8));
//    
//    std::vector<uint32_t> unpackedInput(2*(cbinput/4));
//    std::vector<uint32_t> unpackedOutput(2*(cbinput/4));
//    
//    uint64_t* inputWriteptr = &input[0];
//    uint64_t* inputReadptr = &input[cbinput/8];
//    
//    uint32_t* unpackedInputReadptr = &unpackedInput[cbinput/4];
//    uint32_t* unpackedInputWriteptr = &unpackedInput[0];
//
//    uint32_t* unpackedOutputReadptr = &unpackedOutput[cbinput/4];
//    uint32_t* unpackedOutputWriteptr = &unpackedOutput[0];
//    
//    uint64_t* outputWriteptr = &output[0];
//    uint64_t* outputReadptr = &output[cbinput/8];
//    
//    std::vector<uint32_t> gpuReadOffsets(2*levels+1,0);
//    std::vector<uint32_t> gpuWriteOffsets(2*levels+1,0);
//    
//    std::vector<uint32_t> aboveOverrides(2*levels,0);
//    std::vector<uint32_t> belowOverrides(2*levels,0);
//    
//    int j=0;
//    
//    for (int i=0; i<2*levels+1; i++)
//    {
//        gpuWriteOffsets[i] = j;
//        
//        j+= 2;
//        j = j - (4 & -(j >= 4));
//        j = j + (4 & -(j < 0));
//        
//        gpuReadOffsets[i] = j;
//    }
//    
//    auto cl_instance = init_cl(levels,w,h,bits,"pipeline_kernels.cl",deviceNumber);
//    
//    tbb::task_group group;
//    
//    bool finished = false;
//    
//    int fullness = 0;
//    
//    bool full = false;
//    
//    int tailEnd = 4*levels+5;

//    int image_line = 0;
//    
//    while(1){
//        
//        for (int i=0; i<levels; i++) {
//            if (image_line == 4 + 2*i) aboveOverrides[i] = 0x0 ^ -(levels < 0);
//            else aboveOverrides[i] = 0xFFFFFFFF ^ -(levels < 0);
//            
//            if (image_line == 3 + 2*i) belowOverrides[i] = 0x0 ^ -(levels < 0);
//            else belowOverrides[i] = 0xFFFFFFFF ^ -(levels < 0);
//        }
//        for (int i=levels; i<2*levels; i++) {
//            if (image_line == 4 + 2*i) aboveOverrides[i] = 0xFFFFFFFF ^ -(levels < 0);
//            else aboveOverrides[i] = 0x0 ^ -(levels < 0);
//            
//            if (image_line == 3 + 2*i) belowOverrides[i] = 0xFFFFFFFF ^ -(levels < 0);
//            else belowOverrides[i] = 0x0 ^ -(levels < 0);
//        }
//        
//        group.run([&](){
//            if(!finished && !read_blob(STDIN_FILENO, cbinput, inputWriteptr)) finished = true;
//            
//            unpack_blob_32(cbinput, inputReadptr, unpackedInputWriteptr);
//            
//            pack_blob_32(cbinput, unpackedOutputReadptr, outputWriteptr);
//            
//            if (fullness >= 4*levels+6 || full)
//            {
//                full = true;
//                write_blob(STDOUT_FILENO, cbinput, outputReadptr);
//            }
//            else
//            {
//                fullness++;
//            }
//        });
//        
//        group.run([&](){
//            process_opencl_packed_line(levels, w, bits, gpuReadOffsets, gpuWriteOffsets, unpackedInputReadptr, unpackedOutputWriteptr, aboveOverrides, belowOverrides, cl_instance);
//            
//            for (int i=0; i<2*levels+1; i++)
//            {
//                int j = gpuWriteOffsets[i];
//                
//                j++;
//                j = j - (4 & -(j >= 4));
//                
//                gpuWriteOffsets[i] = j;
//                
//                j = gpuReadOffsets[i];
//                
//                j++;
//                j = j - (4 & -(j >= 4));
//                
//                gpuReadOffsets[i] = j;
//            }
//        });
//        
//        group.wait();
//        
//        if (tailEnd == 0) {
//            break;
//        }
//        if (finished) tailEnd--;
//        
//        std::swap(inputReadptr, inputWriteptr);
//        std::swap(unpackedInputReadptr, unpackedInputWriteptr);
//        std::swap(unpackedOutputReadptr, unpackedOutputWriteptr);
//        std::swap(outputReadptr, outputWriteptr);
//        
//        image_line++;
//        if (image_line == h) image_line = 0;
//    }
//}