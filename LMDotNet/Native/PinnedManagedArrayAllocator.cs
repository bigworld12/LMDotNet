using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet.Native
{
    // from http://blog.getpaint.net/2012/04/30/marshaling-native-arrays-back-as-managed-arrays-without-copying/
    // Note: this allocator is NOT thread safe!
    // TODO: use ConcurrentDictionary (or locks)
    sealed class PinnedManagedArrayAllocator<T> : IDisposable
          where T : struct
    {
        private Dictionary<IntPtr, T[]> baseAddressToManagedArray;
        private Dictionary<IntPtr, GCHandle> baseAddressToGCHandle;
        
        public T[] this[IntPtr pBase] {
            get { return baseAddressToManagedArray[pBase]; }
        }    

        public PinnedManagedArrayAllocator() {
            this.baseAddressToManagedArray = new Dictionary<IntPtr, T[]>();
            this.baseAddressToGCHandle = new Dictionary<IntPtr, GCHandle>();
        }

        // Pass a delegate to this method for DoubleArrayAllocatorDelegate. 
        // Don’t forget to use GC.KeepAlive() on the delegate!
        public IntPtr AllocateArray(int count) {
            T[] managedArray = new T[count];
            var gcHandle = GCHandle.Alloc(managedArray, GCHandleType.Pinned);
            var pBase = gcHandle.AddrOfPinnedObject();
            baseAddressToManagedArray[pBase] = managedArray;
            baseAddressToGCHandle[pBase] = gcHandle;
            return pBase;
        }

        // TODO: implement a correct Dispose() method...
        public void Dispose() {
            if (baseAddressToGCHandle != null) {
                foreach (GCHandle gcHandle in baseAddressToGCHandle.Values) {
                    gcHandle.Free();                    
                }
                baseAddressToGCHandle = null;
            }
            baseAddressToManagedArray = null;            
            GC.SuppressFinalize(this);
        }
    } 
}
