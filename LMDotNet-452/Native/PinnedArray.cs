using System;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    /// <summary>
    /// An array pinned to the heap
    /// </summary>
    /// <typeparam name="T">Type of the elements of the array</typeparam>
    struct PinnedArray<T> where T : struct //:IDisposable
    {
        public readonly T[] Array;
        public readonly GCHandle Handle;
        public readonly IntPtr Address;

        public PinnedArray(T[] array) {
            this.Array = array;
            this.Handle = GCHandle.Alloc(array, GCHandleType.Pinned);
            this.Address = Handle.AddrOfPinnedObject();
        }

        public void Unpin() {
            Handle.Free();
        }

        public static implicit operator T[] (PinnedArray<T> pinnedArray) {
            return pinnedArray.Array;
        }
    }
}
