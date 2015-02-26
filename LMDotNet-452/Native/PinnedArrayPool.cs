using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    // from http://blog.getpaint.net/2012/04/30/marshaling-native-arrays-back-as-managed-arrays-without-copying/
    
    /// <summary>
    /// Allocate arrays and pin them on the managed heap so
    /// they can be used from both managed and native code without
    /// requiring marshalling and/or copying.
    /// When disposed, the pool automatically unpins all allocated arrays.
    /// </summary>
    /// <remarks>This allocator is NOT thread safe!</remarks>    
    /// <typeparam name="T">Type of the array elements; must be an unmanged 
    /// (blitable) type</typeparam>
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

        /// <summary>
        /// Allocate a new pinned T[] and return its base address
        /// </summary>
        /// <param name="count">Length (number of elements) of the new array</param>
        /// <returns>Base address of the pinned array</returns>
        public IntPtr AllocatePinnedArray(int count) {
            var array = new T[count];
            var gch = GCHandle.Alloc(array, GCHandleType.Pinned);
            var pBase = gch.AddrOfPinnedObject();
            pinnedArrays[pBase] = array;
            gcHandles[pBase] = gch;
            return pBase;
        }

        /// <summary>
        /// Unpin the pinned array at the given address
        /// </summary>
        /// <param name="baseAddress">Base address of the array to unpin</param>
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
