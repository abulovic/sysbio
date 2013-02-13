#!/usr/bin/env python

from multiprocessing import Process, Queue

def print_name(name, q):
    print "process ", name
    q.put(name)

def main():
    q = Queue()
    results = []
    for i in range(10):
        p = Process(target=print_name, args=(i,q))
        p.start()
        p.join()

    for k in range(q.qsize()):
        print q.get()
        
if __name__ == '__main__':
    print "before main"
    main()
