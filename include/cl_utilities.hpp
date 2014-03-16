#ifndef HPCE_CW6_GH_cl_utilities_hpp
#define HPCE_CW6_GH_cl_utilities_hpp

#include <vector>
#include <fstream>

#include "CL/cl.hpp"

struct cl_instance
{
    cl::Kernel kernelInstance;
    cl::Buffer proofBuffer;
    cl::CommandQueue commandQueue;
    cl::NDRange offset;
    cl::NDRange globalSize;
    cl::NDRange localSize;
};

std::string LoadSource(const char *fileName)
{
    std::string baseDir="src/kernels";
    if(getenv("HPCE_CL_SRC_DIR")){
        baseDir=getenv("HPCE_CL_SRC_DIR");
    }
    
    std::string fullName=baseDir+"/"+fileName;
    
    std::ifstream src(fullName.c_str(), std::ios::in | std::ios::binary);
    if(!src.is_open())
        throw std::runtime_error("LoadSource : Couldn't load cl file from '"+fullName+"'.");
    
    return std::string((std::istreambuf_iterator<char>(src)),std::istreambuf_iterator<char>());
}

int enumerate_cl_devices()
{
    std::vector<cl::Platform> platforms;
    
	cl::Platform::get(&platforms);
	
    if(platforms.size()==0) return 0;
    
    int selectedPlatform=0;
	if(getenv("HPCE_SELECT_PLATFORM")){
		selectedPlatform=atoi(getenv("HPCE_SELECT_PLATFORM"));
	}
    
	cl::Platform platform=platforms.at(selectedPlatform);
    
    std::vector<cl::Device> devices;
	platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
	
    return devices.size();
}

cl_instance create_cl_instance(std::string source, int deviceNumber = -1)
{
    std::vector<cl::Platform> platforms;
    
	cl::Platform::get(&platforms);
	
    if(platforms.size()==0) throw std::runtime_error("No OpenCL platforms found.");
    
    std::cerr<<"Found "<<platforms.size()<<" platforms\n";
	for(unsigned i=0;i<platforms.size();i++){
		std::string vendor=platforms[0].getInfo<CL_PLATFORM_VENDOR>();
		std::cerr<<"  Platform "<<i<<" : "<<vendor<<"\n";
	}
    
    int selectedPlatform=0;
	if(getenv("HPCE_SELECT_PLATFORM")){
		selectedPlatform=atoi(getenv("HPCE_SELECT_PLATFORM"));
	}
	std::cerr<<"Choosing platform "<<selectedPlatform<<"\n";
	cl::Platform platform=platforms.at(selectedPlatform);
    
    std::vector<cl::Device> devices;
	platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
	if(devices.size()==0){
		throw std::runtime_error("No opencl devices found.\n");
	}
    
	std::cerr<<"Found "<<devices.size()<<" devices\n";
	for(unsigned i=0;i<devices.size();i++){
		std::string name=devices[i].getInfo<CL_DEVICE_NAME>();
		std::cerr<<"  Device "<<i<<" : "<<name<<"\n";
	}
    
    int selectedDevice=0;
    
    if (deviceNumber != -1) selectedDevice = deviceNumber;
    
	std::cerr<<"Choosing device "<<selectedDevice<<"\n";
	cl::Device device=devices.at(selectedDevice);
    
    cl::Context context(devices);
    
    std::string kernelSource=LoadSource(source.c_str());
	
	cl::Program::Sources sources;	// A vector of (data,length) pairs
	sources.push_back(std::make_pair(kernelSource.c_str(), kernelSource.size()+1));	// push on our single string
    
	cl::Program program(context, sources);
	try{
		program.build(devices);
	}catch(...){
		for(unsigned i=0;i<devices.size();i++){
			std::cerr<<"Log for device "<<devices[i].getInfo<CL_DEVICE_NAME>()<<":\n\n";
			std::cerr<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[i])<<"\n\n";
		}
		throw;
	}
    
    //TODO: Determine optimal worker count dynamically.
    // See http://stackoverflow.com/questions/19865127/opencl-ndrange-global-size-local-size
    // and http://stackoverflow.com/questions/3957125/questions-about-global-and-local-work-size
    uint32_t workerCount = 1000;
    
    size_t cbBuffer= 32*workerCount;
    
    cl::Buffer buffer(context, CL_MEM_READ_WRITE, cbBuffer);

    cl::Kernel kernel(program, "bitecoin_miner");
    
    cl::CommandQueue queue(context, device);
    
    cl::NDRange offset(0);                   // Always start iterations at x=0
    cl::NDRange globalSize(workerCount);            // Global size must match the original loops
    cl::NDRange localSize=cl::NullRange;	 // We don't care about local size
    
    return *new cl_instance{kernel,buffer,queue,offset,globalSize,localSize};
}

int test_cl_devices(int levels, unsigned w, unsigned h, unsigned bits, std::string source)
{
    std::vector<cl::Platform> platforms;
    
	cl::Platform::get(&platforms);
	
    if(platforms.size()==0) throw std::runtime_error("No OpenCL platforms found.");
    
    int selectedPlatform=0;
	if(getenv("HPCE_SELECT_PLATFORM")){
		selectedPlatform=atoi(getenv("HPCE_SELECT_PLATFORM"));
	}
    
	cl::Platform platform=platforms.at(selectedPlatform);
    
    std::vector<cl::Device> devices;
	platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
	if(devices.size()==0){
		throw std::runtime_error("No opencl devices found.\n");
	}
    
    uint64_t best = 0xFFFFFFFFFFFFFFFFul;
    int best_device = -1;
    
	for(unsigned i=0;i<devices.size();i++){
        uint64_t microseconds = 0xFFFFFFFFFFFFFFFFul;
        bool failedAttempt = false;
        
        try {
            uint64_t cbinput=uint64_t(w)*bits/8;
            
            std::vector<uint32_t> gpuReadOffsets(2*levels+1,0);
            std::vector<uint32_t> gpuWriteOffsets(2*levels+1,0);
            
            std::vector<uint32_t> aboveOverrides(2*levels,0);
            std::vector<uint32_t> belowOverrides(2*levels,0);
            
            std::vector<uint32_t> unpackedInput(2*(cbinput/4));
            std::vector<uint32_t> unpackedOutput(2*(cbinput/4));
            uint32_t* unpackedInputReadptr = &unpackedInput[cbinput/4];
            uint32_t* unpackedOutputWriteptr = &unpackedOutput[0];
            
            auto cl_instance = create_cl_instance("mining_kernel.cl",i);
            
            timeval before, after;
            
            //gettimeofday(&before, NULL);
            
            //TODO: Add bitecoin timing code here.
            
            for(int i=0; i<100; i++)
            {
                //TODO: Add function to test here
                
                //process_opencl_packed_line(levels, w, bits, gpuReadOffsets, gpuWriteOffsets, unpackedInputReadptr, unpackedOutputWriteptr, aboveOverrides, belowOverrides, cl_instance);
            }
            
            //gettimeofday(&after, NULL);
            //TODO: Add bitecoin timing code here.
            
            
            microseconds =(after.tv_sec - before.tv_sec)*1000000L + after.tv_usec - before.tv_usec;
            
        } catch (std::exception) {
            failedAttempt = true;
        }
        if (failedAttempt) continue;
        if (microseconds < best) {
            best_device = i;
            best = microseconds;
        }
	}
    
    return best_device;
}

#endif
