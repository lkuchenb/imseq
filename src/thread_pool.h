// ============================================================================
// IMSEQ - An immunogenetic sequence analysis tool
// (C) Charite, Universitaetsmedizin Berlin
// Author: Leon Kuchenbecker
// ============================================================================
// 
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published by
// the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
// 
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc., 51
// Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// ============================================================================

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include "thread_check.h"

#ifdef __WITHCDR3THREADS__

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool;

// our worker thread objects
class Worker {
public:
    Worker(ThreadPool &s) : pool(s) { }
    void operator()();
private:
    ThreadPool &pool;
};

// the actual thread pool
class ThreadPool {
public:
    ThreadPool(size_t);
    template<class T, class F>
    std::future<T> enqueue(F f);
    unsigned nWorkers();
    ~ThreadPool();
private:
    friend class Worker;
    
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;
    // the task queue
    std::queue< std::function<void()> > tasks;
    
    // synchronization
    std::mutex queue_mutex, wait_mutex;
    std::condition_variable condition, wait_condition;
    bool stop;
};

// add new work item to the pool
template<class T, class F>
std::future<T> ThreadPool::enqueue(F f)
{
    // don't allow enqueueing after stopping the pool
    if(stop)
        throw std::runtime_error("ThreadPool::enqueue()[E001]");

    auto task = std::make_shared< std::packaged_task<T()> >(f);
    std::future<T> res = task->get_future();    
    { ////////////////////////////////////////////////////////////
        std::unique_lock<std::mutex> lock(queue_mutex);         // 
        tasks.push([task](){ (*task)(); });                     // queue mutex locked
    } ////////////////////////////////////////////////////////////
    condition.notify_one();
    return res;
}

#endif // Multi-threading enabled

#endif // Include guard
