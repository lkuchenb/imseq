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

#include "thread_pool.h"

#ifdef __WITHCDR3THREADS__

unsigned ThreadPool::nWorkers() {
    return workers.size();
}
 
void Worker::operator()()
{
    while(true)
    {
        std::unique_lock<std::mutex> lock(pool.queue_mutex); /////
        while(!pool.stop && pool.tasks.empty())                 //
            pool.condition.wait(lock);                          //
        if(pool.stop && pool.tasks.empty())                     // queue_mutex locked
            return;                                             //
        std::function<void()> task(pool.tasks.front());         //
        pool.tasks.pop();                                       //
        lock.unlock(); ///////////////////////////////////////////
        task();
    }
}
 
// the constructor just launches some amount of workers
ThreadPool::ThreadPool(size_t threads)
    :   stop(false) 
{
    for(size_t i = 0;i<threads;++i)
        workers.push_back(std::thread(Worker(*this)));
}


 
// the destructor joins all threads
ThreadPool::~ThreadPool()
{
    stop = true;
    condition.notify_all();
    for(size_t i = 0;i<workers.size();++i)
        workers[i].join();
}

#endif 
