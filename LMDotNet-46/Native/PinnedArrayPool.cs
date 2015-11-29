using System;
using System.Collections.Generic;
using System.Diagnostics;

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
        private Dictionary<IntPtr, PinnedArray<T>> pool;
        
        public T[] this[IntPtr pBase] {
            get { return pool[pBase].Array; }
        }    

        public PinnedArrayPool() {
            this.pool = new Dictionary<IntPtr, PinnedArray<T>>();            
        }

        /// <summary>
        /// Allocate a new pinned T[] and return a PinnedArray instance;
        /// supposed to be callaed from managed code
        /// </summary>
        /// <param name="count">Length (number of elements) of the new array</param>
        /// <returns>Base address of the pinned array</returns>
        public PinnedArray<T> AllocatePinnedArray(int count) {
            var pinnedArray = new PinnedArray<T>(new T[count]);
            pool[pinnedArray.Address] = pinnedArray;
            return pinnedArray;
        }

        /// <summary>
        /// Allocate a new pinned T[] and return its base address,
        /// supposed to be callaed from unmanaged code as an allocator
        /// ("calloc")
        /// </summary>
        /// <param name="count">Length (number of elements) of the new array</param>
        /// <returns>Base address of the pinned array</returns>
        public IntPtr Calloc(int count) {
            var pinnedArray = AllocatePinnedArray(count);
            return pinnedArray.Address;
        }

        /// <summary>
        /// Unpin the given pinned array and remove it
        /// from the pool (deallocate)
        /// </summary>
        /// <param name="pinnedArray">The pinned array to unpin and remove from the pool</param>
        public void UnpinArray(PinnedArray<T> pinnedArray) {
            pinnedArray.Unpin();
            pool.Remove(pinnedArray.Address);
        }

        /// <summary>
        /// Unpin the pinned array at the given address,
        /// supposed to be callaed from unmanaged code as an deallocator
        /// ("free")
        /// </summary>
        /// <param name="baseAddress">Base address of the array to free</param>
        public void Free(IntPtr baseAddress) {
            UnpinArray(pool[baseAddress]);
        }

        // TODO: implement a correct Dispose() method...
        public void Dispose() {
            if (pool != null) {
                foreach (var pinnedArr in pool.Values) {
                    pinnedArr.Unpin();
                }
                pool = null;
            }
            pool = null;            
            GC.SuppressFinalize(this);
        }
    } 
}
