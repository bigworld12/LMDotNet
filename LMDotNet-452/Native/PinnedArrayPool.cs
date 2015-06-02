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
        struct PinnedArray
        {
            public readonly T[] Array;
            public readonly GCHandle Handle;

            public PinnedArray(T[] array, GCHandle handle) {
                this.Array = array;
                this.Handle = handle;
            }

            //public static Pin
            //public Free
            //public implicit
        }

        private Dictionary<IntPtr, PinnedArray> pool;
        
        public T[] this[IntPtr pBase] {
            get { return pool[pBase].Array; }
        }    

        public PinnedArrayPool() {
            this.pool = new Dictionary<IntPtr, PinnedArray>();            
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
            pool[pBase] = new PinnedArray(array, gch);
            return pBase;
        }

        /// <summary>
        /// Unpin the pinned array at the given address
        /// </summary>
        /// <param name="baseAddress">Base address of the array to unpin</param>
        public void UnpinArray(IntPtr baseAddress) {
            pool[baseAddress].Handle.Free();
            pool.Remove(baseAddress);
        }

        // TODO: implement a correct Dispose() method...
        public void Dispose() {
            if (pool != null) {
                foreach (var arr in pool.Values) {
                    arr.Handle.Free();
                }
                pool = null;
            }
            pool = null;            
            GC.SuppressFinalize(this);
        }
    } 
}
