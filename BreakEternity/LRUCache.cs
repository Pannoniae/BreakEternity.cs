using System;
using System.Collections.Generic;

namespace BreakEternity {

    public class LRUCache<K, V> where K : notnull {
        private readonly Dictionary<K, ListNode<K, V>?> map = new();
        private ListNode<K, V>? first;
        private ListNode<K, V>? last;
        public readonly double maxSize;

        /**
   * @param maxSize The maximum size for this cache. We recommend setting this
   * to be one less than a power of 2, as most hashtables - including V8's
   * Object hashtable (https://crsrc.org/c/v8/src/objects/ordered-hash-table.cc)
   * - uses powers of two for hashtable sizes. It can't exactly be a power of
   * two, as a .set() call could temporarily set the size of the map to be
   * maxSize + 1.
   */
        public LRUCache(double maxSize) {
            this.maxSize = maxSize;
        }

        public double Count => map.Count;


        /**
   * Gets the specified key from the cache, or undefined if it is not in the
   * cache.
   * @param key The key to get.
   * @returns The cached value, or undefined if key is not in the cache.
   */
        public V? get(K key) {
            var node = map[key];
            if (node == null) {
                return default;
            }
            // It is guaranteed that there is at least one item in the cache.
            // Therefore, first and last are guaranteed to be a ListNode...
            // but if there is only one item, they might be the same.

            // Update the order of the list to make this node the first node in the
            // list.
            // This isn't needed if this node is already the first node in the list.
            if (node != first) {
                // As this node is DIFFERENT from the first node, it is guaranteed that
                // there are at least two items in the cache.
                // However, this node could possibly be the last item.
                if (node == last) {
                    // This node IS the last node.
                    last = node.prev;
                    // From the invariants, there must be at least two items in the cache,
                    // so node - which is the original "last node" - must have a defined
                    // previous node. Therefore, this.last - set above - must be defined
                    // here.
                    last!.next = null;
                }
                else {
                    // This node is somewhere in the middle of the list, so there must be at
                    // least THREE items in the list, and this node's prev and next must be
                    // defined here.
                    node.prev!.next = node.next;
                    node.next!.prev = node.prev;
                }

                node.next = first;
                // From the invariants, there must be at least two items in the cache, so
                // this.first must be a valid ListNode.
                first!.prev = node;
                first = node;
            }

            return node.value;
        }

        /**
   * Sets an entry in the cache.
   *
   * @param key The key of the entry.
   * @param value The value of the entry.
   * @throws Error, if the map already contains the key.
   */
        public void set(K key, V value) {
            // Ensure that this.maxSize >= 1.
            if (maxSize < 1) {
                return;
            }

            if (map.ContainsKey(key)) {
                throw new ArgumentException("Cannot update existing keys in the cache");
            }

            var node = new ListNode<K, V>(key, value);
            // Move node to the front of the list.
            if (first == null) {
                // If the first is undefined, the last is undefined too.
                // Therefore, this cache has no items in it.
                first = node;
                last = node;
            }
            else {
                // This cache has at least one item in it.
                node.next = first;
                first.prev = node;
                first = node;
            }

            map[key] = node;

            while (map.Count > maxSize) {
                // We are guaranteed that this.maxSize >= 1,
                // so this.map.size is guaranteed to be >= 2,
                // so this.first and this.last must be different valid ListNodes,
                // and this.last.prev must also be a valid ListNode (possibly this.first).
                var lastNode = this.last!;
                map.Remove(lastNode.key);
                last = lastNode.prev;
                last!.next = null;
            }
        }
    }

    /**
 * A node in a doubly linked list.
 */
    public class ListNode<K, V> {
        public K key;
        public V value;
        public ListNode<K, V>? next;
        public ListNode<K, V>? prev;

        public ListNode(K key, V value) {
            this.key = key;
            this.value = value;
        }
    }
}