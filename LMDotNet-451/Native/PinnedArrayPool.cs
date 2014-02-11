using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    // from http://blog.getpaint.net/2012/04/30/marshaling-native-arrays-back-as-managed-arrays-without-copying/
    // Note: this allocator is NOT thread safe!
    [DebuggerDisplay("array count = {pinnedArrays.Count}")]
    sealed class PinnedArrayPool<T> : IDisposable
          where T : struct
    {
        private Dictionary<IntPtr, T[]> pinnedArrays;
        private Dictionary<IntPtr, GCHandle> gcHandles;
        
        public T[] this[IntPtr pBase] {
            get { return pinnedArrays[pBase]; }
        }    

        public PinnedArrayPool() {
            this.pinnedArrays = new Dictionary<IntPtr, T[]>();
            this.gcHandles = new Dictionary<IntPtr, GCHandle>();
        }

        public IntPtr AllocatePinnedArray(int count) {
            var array = new T[count];
            var gch = GCHandle.Alloc(array, GCHandleType.Pinned);
            var pBase = gch.AddrOfPinnedObject();
            pinnedArrays[pBase] = array;
            gcHandles[pBase] = gch;
            return pBase;
        }

        // un-pin array
        public void UnpinArray(IntPtr baseAddress) {
            gcHandles[baseAddress].Free();
            gcHandles.Remove(baseAddress);
            pinnedArrays.Remove(baseAddress);
        }

        // TODO: implement a correct Dispose() method...
        public void Dispose() {
            if (gcHandles != null) {
                foreach (var gch in gcHandles.Values) {
                    gch.Free();                    
                }
                gcHandles = null;
            }
            pinnedArrays = null;            
            GC.SuppressFinalize(this);
        }
    } 
}
