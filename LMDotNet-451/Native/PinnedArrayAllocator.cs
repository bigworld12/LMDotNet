using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    // from http://blog.getpaint.net/2012/04/30/marshaling-native-arrays-back-as-managed-arrays-without-copying/
    // Note: this allocator is NOT thread safe!
    // TODO: use ConcurrentDictionary (or locks)
    [DebuggerDisplay("array count = {ptrToArray.Count}")]
    sealed class PinnedArrayAllocator<T> : IDisposable
          where T : struct
    {
        private Dictionary<IntPtr, T[]> ptrToArray;
        private Dictionary<IntPtr, GCHandle> ptrToHandle;
        
        public T[] this[IntPtr pBase] {
            get { return ptrToArray[pBase]; }
        }    

        public PinnedArrayAllocator() {
            this.ptrToArray = new Dictionary<IntPtr, T[]>();
            this.ptrToHandle = new Dictionary<IntPtr, GCHandle>();
        }

        // Pass a delegate to this method for DoubleArrayAllocatorDelegate. 
        // Don’t forget to use GC.KeepAlive() on the delegate!
        public IntPtr AllocatePinnedArray(int count) {
            var array = new T[count];
            var gch = GCHandle.Alloc(array, GCHandleType.Pinned);
            var pBase = gch.AddrOfPinnedObject();
            ptrToArray[pBase] = array;
            ptrToHandle[pBase] = gch;
            return pBase;
        }

        // un-pin array
        public void UnpinArray(IntPtr baseAddress) {
            ptrToHandle[baseAddress].Free();
            ptrToHandle.Remove(baseAddress);
            ptrToArray.Remove(baseAddress);
        }

        // TODO: implement a correct Dispose() method...
        public void Dispose() {
            if (ptrToHandle != null) {
                foreach (var gch in ptrToHandle.Values) {
                    gch.Free();                    
                }
                ptrToHandle = null;
            }
            ptrToArray = null;            
            GC.SuppressFinalize(this);
        }
    } 
}
