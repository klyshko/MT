#pragma once

#include <iostream>
#include <string>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <vector_functions.h>
#include <cuda_runtime_api.h>
#include "wrapper.h"

#ifdef __CUDACC__
#define CUDA_FUNC __device__ __host__
#define cuda_assert(x)
#else
#define CUDA_FUNC
#define cuda_assert assert
#endif

#ifndef HALFWARP_SIZE
#define HALFWARP_SIZE 16
#endif

#define CUDA_SIMPLE_INVOKE_BLOCK_SIZE (4*8*HALFWARP_SIZE)

#define CUDA_STATIC_ASSERT(cond) typedef char _StaticAssertFailed_[(cond)? 1: -1]

// NOTE: this implementation is not safe, bacause it consist of two expressions
#define CUDA_INVOKE_BA(kernel, block_dim, thread_dim, before_code, after_code) \
struct _CUDA_INVOKE_##kernel { \
							   _CUDA_INVOKE_##kernel() { before_code } \
							   ~_CUDA_INVOKE_##kernel() { after_code } \
						   }; \
kernel<<<(_CUDA_INVOKE_##kernel(), cudaCheckBlockDim(block_dim, #kernel, __FILE__, __LINE__)), cudaCheckThreadDim(thread_dim, #kernel, __FILE__, __LINE__)>>> \
/**/

#ifdef CUDA_STRICT
#define CUDA_INVOKE(kernel, block_dim, thread_dim) \
kernel<<<(_StrictInvokeChecker_(#kernel, __FILE__, __LINE__), cudaCheckBlockDim(block_dim, #kernel, __FILE__, __LINE__)), \
cudaCheckThreadDim(thread_dim, #kernel, __FILE__, __LINE__)>>>
#else
#define CUDA_INVOKE(kernel, block_dim, thread_dim) \
kernel<<<cudaCheckBlockDim(block_dim, #kernel, __FILE__, __LINE__), \
cudaCheckThreadDim(thread_dim, #kernel, __FILE__, __LINE__)>>>
#endif

#define CUDA_SIMPLE_INVOKE(kernel, total_threads) CUDA_INVOKE(kernel, (unsigned)ceil((float)(total_threads)/CUDA_SIMPLE_INVOKE_BLOCK_SIZE), CUDA_SIMPLE_INVOKE_BLOCK_SIZE)

#define CUDA_INVOKE_THREADS(kernel, total_threads, threads_per_block) CUDA_INVOKE(kernel, (unsigned)ceil((float)(total_threads)/threads_per_block), threads_per_block)

class CudaException: std::runtime_error {
public:
	explicit CudaException(cudaError_t error, const char* kernel_name = 0, const char* file = 0, int line = 0)
		: std::runtime_error(FormatWhat(error, kernel_name, file, line).c_str())
		, error_(error)
		, kernel_name_(kernel_name)
		, file_(file)
		, line_(line)
	{}
	~CudaException() throw() {}
	cudaError_t error() const { return error_; }
	static std::string FormatWhat(cudaError_t error, const char* kernel_name, const char* file, int line) {
		char buf[1024];
		std::string ret =  + "CudaException: " + std::string(cudaGetErrorString(error));
		if (kernel_name)
			ret = "[" + std::string(kernel_name) + "]:" + ret;
		if (line) {
			sprintf(buf, "%d", line);
			ret = std::string(buf) + ":" + ret;
		}
		if (file)
			ret = std::string(file) + ":" + ret;
		return ret;
	}

private:
	cudaError_t error_;
	const char* kernel_name_;
	const char* file_;
	int line_;
};

// WARNING: This function is not thread-safe!!!
//          Use it only if you have the only thread working with CUDA
inline const cudaDeviceProp& cudaCurrentDeviceProp() {
	static bool first_call = true;
	static int device;
	static cudaDeviceProp prop;
	if (first_call) {
		cudaError_t error = cudaGetDevice(&device);
		if (error != cudaSuccess)
			throw CudaException(error);
		error = cudaGetDeviceProperties(&prop, device);
		if (error != cudaSuccess)
			throw CudaException(error);
		first_call = false;
	}
	return prop;
}

inline dim3 cudaCheckBlockDim(const dim3& dim, const char* kernel_name, const char* file, int line) {
	const cudaDeviceProp& prop = cudaCurrentDeviceProp();
	if ((dim.x > (unsigned)prop.maxGridSize[0] || dim.y > (unsigned)prop.maxGridSize[1] || dim.z > (unsigned)prop.maxGridSize[2])||
		(dim.x == 0) || (dim.y == 0) || (dim.z == 0))
		throw CudaException(cudaErrorLaunchFailure, kernel_name, file, line);
	return dim;
}

inline dim3 cudaCheckThreadDim(const dim3& dim, const char* kernel_name, const char* file, int line) {
	const cudaDeviceProp& prop = cudaCurrentDeviceProp();
	if (dim.x > (unsigned)prop.maxThreadsDim[0] || dim.y > (unsigned)prop.maxThreadsDim[1] || dim.z > (unsigned)prop.maxThreadsDim[2] ||
		dim.x * dim.y * dim.z > (unsigned)prop.maxThreadsPerBlock || (dim.x == 0) || (dim.y == 0) || (dim.z == 0))
		throw CudaException(cudaErrorLaunchFailure, kernel_name, file, line);
	return dim;
}

#define cudaProcessError() cudaProcessErrorImpl(__FILE__, __LINE__, 0)
#define checkCUDAError(msg) cudaProcessErrorImpl(__FILE__, __LINE__, 0, msg)

inline void cudaProcessErrorImpl(const char* file, int line, const char* kernel_name, const char* msg = "") {
	//cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error == cudaSuccess)
		error = cudaThreadSynchronize();
	if (error == cudaSuccess)
		error = cudaGetLastError();
	if (error != cudaSuccess) {
		if (kernel_name)
			die("%s:%d [%s]: CudaError: %s: %s\n", file, line, kernel_name, msg, cudaGetErrorString(error));
		else
			die("%s:%d: CudaError: %s: %s\n", file, line, msg, cudaGetErrorString(error));
	}
}

void cudaSetDefaultDevice();

struct _StrictInvokeChecker_ {
	const char* kernel_name;
	const char* file;
	int line;
	_StrictInvokeChecker_(const char* k, const char* f, int l)
		: kernel_name(k), file(f), line(l) {
		cudaProcessErrorImpl(file, line, ("before " + std::string(kernel_name)).c_str());
	}
	~_StrictInvokeChecker_() {
		cudaProcessErrorImpl(file, line, ("after " + std::string(kernel_name)).c_str());
	}
};

//! Image of data stored on device. A copy is created on stack.
//! Syntax for working with HostImage objects is identical to ordinary pointers syntax.
template <class T>
class HostImage {
public:
	explicit HostImage(T* d_data)
		: d_data_(d_data) {
		Load();
	}

	void Load() {
		cudaMemcpy(&image_, d_data_, sizeof(T), cudaMemcpyDeviceToHost);
	}

	void Save() {
		cudaMemcpy(d_data_, &image_, sizeof(T), cudaMemcpyHostToDevice);
	}

	// Image access
	T* operator->() { return &image_; }
	const T* operator->() const { return &image_; }
	T& operator*() { return image_; }
	const T& operator*() const { return image_; }
	operator T*() { return &image_; }
	operator const T*() const { return &image_; }

	// Device data pointer access
	T* d_data() { return d_data_; }
	const T* d_data() const { return d_data_; }

private:
	T image_;
	T* d_data_;
};

template <class T>
inline HostImage<T> ToHost(T* data) {
	return HostImage<T>(data);
}


CUDA_FUNC inline float3 make_float3(const float4 &v) {
	return make_float3(v.x, v.y, v.z);
}

// Helpers for built-in cuda vector types
CUDA_FUNC inline float findDistance(float4 a, float4 b) {
	return sqrtf((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

CUDA_FUNC inline float4 operator+(const float4& a, const float4& b) {
	return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

CUDA_FUNC inline float2 operator+(const float2& a, const float2& b) {
	return make_float2(a.x + b.x, a.y + b.y);
}

CUDA_FUNC inline float4 operator-(const float4& a, const float4& b) {
	return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

CUDA_FUNC inline float2 operator-(const float2& a, const float2& b) {
	return make_float2(a.x - b.x, a.y - b.y);
}

CUDA_FUNC inline void operator+=(float4& a, const float4& b) {
	a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
}

CUDA_FUNC inline void operator+=(float2& a, const float2& b) {
	a.x += b.x; a.y += b.y;
}

CUDA_FUNC inline void operator-=(float4& a, const float4& b) {
	a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
}

CUDA_FUNC inline void operator-=(float2& a, const float2& b) {
	a.x -= b.x; a.y -= b.y;
}

CUDA_FUNC inline void operator/=(float4& a, float b) {
	a.x /= b; a.y /= b; a.z /= b; a.w /= b;
}

CUDA_FUNC inline void operator/=(float2& a, float b) {
	a.x /= b; a.y /= b;
}

CUDA_FUNC inline float4 operator/(const float4& a, float b) {
	return make_float4(a.x / b, a.y / b, a.z / b, a.w / b);
}

CUDA_FUNC inline float2 operator/(const float2& a, float b) {
	return make_float2(a.x / b, a.y / b);
}

CUDA_FUNC inline void operator*=(float4& a, float b) {
	a.x *= b; a.y *= b; a.z *= b; a.w *= b;
}

CUDA_FUNC inline void operator*=(float2& a, float b) {
	a.x *= b; a.y *= b;
}

CUDA_FUNC inline float4 operator*(const float4& a, float b) {
	return make_float4(a.x * b, a.y * b, a.z * b, a.w * b);
}

CUDA_FUNC inline float2 operator*(const float2& a, float b) {
	return make_float2(a.x * b, a.y * b);
}

CUDA_FUNC inline float dot(const float4 &a, const float4 &b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

CUDA_FUNC inline float4 cross(const float4 &a, const float4 &b) {
	return make_float4(a.y*b.z-b.y*a.z, a.z*b.x-b.z*a.x, a.x*b.y-b.x*a.y, 0.0f);
}

CUDA_FUNC inline float abs2(const float4 &a) {
	return a.x * a.x + a.y * a.y + a.z * a.z;
}

CUDA_FUNC inline float abs(const float4 &a) {
	return sqrtf(abs2(a));
}

CUDA_FUNC inline float normalize(float4 &a) {
	float v = abs(a);
	a /= v;
	return v;
}

CUDA_FUNC inline float normalize_s(float4 &a) { // Safe normalize
	float v = abs(a);
	if (v < 1e-5) // Chosen arbitrary. It just seemed suitable
		a = make_float4(1.0f, 0.0f, 0.0f, 0.0f); // Chosen arbitrary
	else
		a /= v;
	return v;
}

// Aargh, `isnan' is actually a macro, so no overloading :-(
// I've read somewhere that it will be made templaate in C++0x, but who cares
// Altogether, this piece of code sucks, but I'm f-cking lazy
CUDA_FUNC inline int my_isnan(const float4 &a) {
	return (isnan(a.x) || isnan(a.y) || isnan(a.z) || isnan(a.w));
}

CUDA_FUNC inline int my_isnan(const float3 &a) {
	return (isnan(a.x) || isnan(a.y) || isnan(a.z));
}

CUDA_FUNC inline int my_isnan(const float2 &a) {
		return (isnan(a.x) || isnan(a.y));
}

CUDA_FUNC inline int my_isnan(const float &a) {
	return isnan(a);
}

/*
CUDA_FUNC float4 make_float4(const float3 &a, float w) {
    return make_float4(a.x,a.y,a.z,w);
}*/

inline std::ostream& operator << (std::ostream& os, int4 v) {
	return os << "[" << v.x << ", " << v.y << ", " << v.z << "; " << v.w << "]";
}

inline std::ostream& operator << (std::ostream& os, float3 v) {
	return os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
}

inline std::ostream& operator << (std::ostream& os, float4 v) {
	return os << "[" << v.x << ", " << v.y << ", " << v.z << "; " << v.w << "]";
}
