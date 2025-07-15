import os
import sys
import pickle
import argparse

class Sequence:
    def __init__(self, path: str) -> None:
        self.path = path

    def read(self) -> str:
        """Read DNA sequence from file, skipping FASTA headers."""
        seq = ""
        with open(self.path, 'r') as file:
            for line in file:
                if not line.startswith('>'):
                    seq += line.strip().upper()
        return seq

    def write(self, content: str) -> None:
        """Write content to file."""
        with open(self.path, 'w') as file:
            file.write(content)

class SuffixArray:
    @staticmethod
    def build_suffix_array(text: str) -> list:
        """
        Build suffix array for the given text using a naive but clear approach.
        
        Args:
            text: The input string
            
        Returns:
            List of indices representing the sorted order of suffixes
        """
        # Generate all suffixes with their original indices
        suffixes = [(i, text[i:]) for i in range(len(text))]
        
        # Sort suffixes lexicographically
        suffixes.sort(key=lambda x: x[1])
        
        # Extract indices to form the suffix array
        return [idx for idx, _ in suffixes]
    
    @staticmethod
    def build_suffix_array_efficient(text: str) -> list:
        """
        Build suffix array using a more efficient algorithm.
        This is a prefix doubling implementation with O(n log n) time complexity.
        
        Args:
            text: The input string
            
        Returns:
            List of indices representing the sorted order of suffixes
        """
        n = len(text)
        
        # Initial ranking of characters
        ranks = [ord(c) for c in text]
        
        # Add a sentinel value for the "empty" suffix
        ranks.append(-1)
        
        # Initialize suffix array
        sa = list(range(n))
        
        # Sort suffixes based on first character
        sa.sort(key=lambda i: ranks[i])
        
        # Array to store new ranks
        new_ranks = [0] * (n + 1)
        new_ranks[n] = -1  # Sentinel rank
        
        # Prefix doubling
        h = 1
        while h < n:
            # Sort by rank pairs (current, current+h)
            sa.sort(key=lambda i: (ranks[i], ranks[min(i + h, n)]))
            
            # Update ranks based on sorted order
            new_ranks[sa[0]] = 0
            for i in range(1, n):
                # Compare current suffix with previous suffix
                prev = sa[i-1]
                curr = sa[i]
                
                # Check if they have the same rank pair
                if (ranks[prev], ranks[min(prev + h, n)]) == (ranks[curr], ranks[min(curr + h, n)]):
                    new_ranks[curr] = new_ranks[prev]  # Same rank as previous
                else:
                    new_ranks[curr] = new_ranks[prev] + 1  # New rank
            
            # Update ranks
            ranks = new_ranks.copy()
            
            # If all suffixes have unique ranks, we're done
            if ranks[sa[n-1]] == n-1:
                break
                
            # Double the prefix length
            h *= 2
            
        return sa

class BurrowsWheeler:
    @staticmethod
    def bwt_with_suffix_array(sequence: str) -> str:
        """
        Burrows-Wheeler Transform implementation using suffix arrays.
        
        Args:
            sequence: The input string to transform
            
        Returns:
            The BWT of the input string
        """
        # Ensure the sequence ends with a terminator
        if not sequence.endswith('$'):
            sequence += '$'
        
        n = len(sequence)
        
        # Build suffix array
        sa = SuffixArray.build_suffix_array_efficient(sequence)
        
        # Compute BWT
        bwt = []
        for i in sa:
            bwt.append(sequence[(i - 1) % n])  # Get the character before each suffix
        
        return ''.join(bwt)

    @staticmethod
    def merge_sort(arr):
        if len(arr) <= 1:
            return arr

        mid = len(arr) // 2
        left = BurrowsWheeler.merge_sort(arr[:mid])
        right = BurrowsWheeler.merge_sort(arr[mid:])

        return BurrowsWheeler.merge(left, right)

    @staticmethod
    def merge(left, right):
        merged = []
        i = j = 0

        while i < len(left) and j < len(right):
            if left[i] < right[j]:  # Lexicographical comparison
                merged.append(left[i])
                i += 1
            else:
                merged.append(right[j])
                j += 1

        merged.extend(left[i:])
        merged.extend(right[j:])
        return merged

    @staticmethod
    def bwt(sequence: str) -> str:
        """
        Burrows-Wheeler Transform using merge sort instead of built-in sorting.
        """
        if not sequence.endswith('$'):
            sequence += '$'
        n = len(sequence)
        rotations = [sequence[i:] + sequence[:i] for i in range(n)]
        sorted_rotations = BurrowsWheeler.merge_sort(rotations)
        return ''.join([rot[-1] for rot in sorted_rotations])


    @staticmethod
    def inverse_bwt(bwt_string: str) -> str: 
        """Inverse Burrows-Wheeler Transform."""
        # Check if BWT string is valid
        if not bwt_string or '$' not in bwt_string:
            print(f"WARNING: Invalid BWT string for inverse: '{bwt_string}'")
            return ""
        
        # Get the length of the string
        n = len(bwt_string)
        
        # Create sorted characters array (first column)
        sorted_chars = sorted([(i, char) for i, char in enumerate(bwt_string)], key=lambda x: x[1])
        
        # Track the original position in the BWT string
        result = []
        next_idx = bwt_string.index('$')  # Start with the row containing $
        
        # Follow the cycle for n-1 steps (exclude the $ at the end)
        for _ in range(n - 1):
            next_idx = sorted_chars[next_idx][0]
            result.append(bwt_string[next_idx])
        
        # Return reconstructed string
        return ''.join(result)

class MoveToFront:
    @staticmethod
    def encode(text: str) -> list:
        """
        Apply Move-to-Front encoding on a string.
        
        Args:
            text: The input string to encode
            
        Returns:
            A list of integer indices representing the encoded sequence
        """
        # Create alphabet with unique characters from text
        alphabet = sorted(set(text))
        result = []
        
        # Process each character in the input text
        for char in text:
            # Find position of character in current alphabet
            index = alphabet.index(char)
            result.append(index)
            
            # Move the character to front
            alphabet.pop(index)
            alphabet.insert(0, char)
            
        return result
    
    @staticmethod
    def decode(indices: list, alphabet: list) -> str:
        """
        Apply Move-to-Front decoding on a list of indices.
        
        Args:
            indices: List of integer indices from MTF encoding
            alphabet: Initial alphabet (sorted unique characters from original text)
            
        Returns:
            The decoded string
        """
        result = []
        
        # Create a copy of the alphabet to work with
        working_alphabet = alphabet.copy()
        
        # Process each index
        for idx in indices:
            # Get character at current index
            char = working_alphabet[idx]
            result.append(char)
            
            # Move the character to front
            working_alphabet.pop(idx)
            working_alphabet.insert(0, char)
            
        return ''.join(result)

class HuffmanNode:
    def __init__(self, char: str = None, freq: int = 0, left=None, right=None):
        self.char = char
        self.freq = freq
        self.left = left
        self.right = right

class HuffmanCoding:
    def __init__(self, sequence: list = None, codes: dict = None):
        self.codes = codes or {}
        if sequence is not None:  # Changed to handle sequence as list of integers
            self.sequence = sequence
            self.root = self.build_tree()
            self.generate_codes(self.root)
        elif codes:
            self.root = self.rebuild_tree_from_codes()
        else:
            raise ValueError("Must provide either sequence or codes")

    def build_tree(self):
        """Build Huffman tree based on character/index frequencies."""
        freq = {}
        for value in self.sequence:
            freq[value] = freq.get(value, 0) + 1

        nodes = [HuffmanNode(char, freq) for char, freq in freq.items()]

        while len(nodes) > 1:
            nodes = sorted(nodes, key=lambda x: x.freq)
            left = nodes.pop(0)
            right = nodes.pop(0)
            merged = HuffmanNode(None, left.freq + right.freq, left, right)
            nodes.append(merged)

        return nodes[0] if nodes else None

    def rebuild_tree_from_codes(self) -> HuffmanNode:
        """Rebuild Huffman tree from existing codes."""
        root = HuffmanNode()
        for char, code in self.codes.items():
            current = root
            for bit in code:
                if bit == '0':
                    if not current.left:
                        current.left = HuffmanNode()
                    current = current.left
                else:
                    if not current.right:
                        current.right = HuffmanNode()
                    current = current.right
            current.char = char
        return root

    def generate_codes(self, node, current_code=""):
        """Generate Huffman codes for each character/index."""
        if node is None:
            return

        if node.char is not None:
            self.codes[node.char] = current_code
            return

        self.generate_codes(node.left, current_code + "0")
        self.generate_codes(node.right, current_code + "1")

    def encode(self) -> str:
        """Encode the sequence using generated Huffman codes."""
        return ''.join([self.codes[value] for value in self.sequence])

    def decode(self, encoded_str: str) -> list:
        """Decode Huffman encoded string back to original MTF indices."""
        current_node = self.root
        decoded = []
        
        for bit in encoded_str:
            if bit == '0':
                current_node = current_node.left
            else:
                current_node = current_node.right

            if current_node is None:
                raise ValueError("Invalid Huffman code encountered")

            if current_node.char is not None:
                decoded.append(current_node.char)
                current_node = self.root

        if current_node != self.root:
            raise ValueError("Incomplete Huffman code at end of input")

        return decoded

class DNACompressor:
    def __init__(self, input_path: str, chunk_size: int = 1000000, use_suffix_array: bool = True):
        self.input_path = input_path
        self.chunk_size = chunk_size
        self.output_dir = os.path.dirname(input_path)
        self.base_name = os.path.splitext(os.path.basename(input_path))[0]
        self.use_suffix_array = use_suffix_array

    def compress(self, output_path: str) -> str:
        """Compress DNA sequence using BWT + MTF + Huffman coding."""
        print(f"Compressing {self.input_path}...")
        
        # Read input DNA sequence
        seq = Sequence(self.input_path).read()
        if not seq:
            raise ValueError("Input file is empty or contains no sequence data")
        
        print(f"Input sequence: '{seq[:20]}{'...' if len(seq) > 20 else ''}'")
        
        # Apply BWT with suffix array optimization
        print(f"Applying Burrows-Wheeler Transform {'with suffix array' if self.use_suffix_array else '(original method)'}...")
        
        if self.use_suffix_array:
            bwt_result = BurrowsWheeler.bwt_with_suffix_array(seq)
        else:
            bwt_result = BurrowsWheeler.bwt(seq)
            
        print(f"BWT result: '{bwt_result[:20]}{'...' if len(bwt_result) > 20 else ''}'")
        
        # Apply Move-to-Front encoding
        print("Applying Move-to-Front encoding...")
        mtf_result = MoveToFront.encode(bwt_result)
        print(f"MTF result (first 20 values): {mtf_result[:20]}{'...' if len(mtf_result) > 20 else ''}")
        
        # Apply Huffman coding
        print("Applying Huffman Coding...")
        huffman = HuffmanCoding(mtf_result)
        encoded_bits = huffman.encode()
        print(f"Huffman codes: {huffman.codes}")
        print(f"Encoded bits (first 100 bits): {encoded_bits[:100]}{'...' if len(encoded_bits) > 100 else ''}")
        
        # Save compressed data and metadata
        metadata = {
            'huffman_codes': huffman.codes,
            'original_length': len(seq),
            'bwt_length': len(bwt_result),
            'alphabet': sorted(set(bwt_result)),  # Save alphabet for MTF decoding
            'use_suffix_array': self.use_suffix_array,
            'bwt_result': bwt_result,  # For debugging, can be removed in production
            'mtf_result': mtf_result   # For debugging, can be removed in production
        }
        
        # Write binary output as text (0s and 1s)
        with open(output_path, 'w') as f:
            f.write(encoded_bits)
        
        # Save metadata
        metadata_path = output_path + '.meta'
        with open(metadata_path, 'wb') as f:
            pickle.dump(metadata, f)
        
        print(f"Compression complete. Output saved to {output_path}")
        print(f"Original size: {len(seq)} bytes")
        compressed_size_bytes = len(encoded_bits) / 8  # Convert bits to bytes
        print(f"Compressed size: {compressed_size_bytes:.2f} bytes (as bits: {len(encoded_bits)})")
        if len(seq) > 0:
            ratio = len(seq) / compressed_size_bytes if compressed_size_bytes > 0 else float('inf')
            print(f"Compression ratio: {ratio:.2f}x")
        return output_path

    def decompress(self, compressed_path: str, output_path: str) -> str:
        """Decompress binary file back to DNA sequence."""
        print(f"Decompressing {compressed_path}...")
        
        # Load metadata
        metadata_path = compressed_path + '.meta'
        try:
            with open(metadata_path, 'rb') as f:
                metadata = pickle.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Metadata file not found: {metadata_path}\n"
                "Please ensure both the compressed file and its .meta file are present"
            )
        
        print(f"Metadata loaded with keys: {list(metadata.keys())}")
        
        # Read compressed binary data as text
        with open(compressed_path, 'r') as f:
            encoded_bits = f.read()
        
        print(f"Encoded bits read (first 100): {encoded_bits[:100]}{'...' if len(encoded_bits) > 100 else ''}")
        
        # Rebuild Huffman tree from codes
        print("Rebuilding Huffman tree...")
        try:
            huffman = HuffmanCoding(codes=metadata['huffman_codes'])
            print(f"Huffman codes from metadata: {list(metadata['huffman_codes'].items())[:5]}{'...' if len(metadata['huffman_codes']) > 5 else ''}")
        except Exception as e:
            raise ValueError(f"Failed to rebuild Huffman tree: {str(e)}")
        
        # Decode Huffman
        print("Decoding Huffman...")
        try:
            mtf_indices = huffman.decode(encoded_bits)
            print(f"Decoded MTF indices (first 20): {mtf_indices[:20]}{'...' if len(mtf_indices) > 20 else ''}")
            
            # Verify MTF indices if available for debugging
            if 'mtf_result' in metadata and mtf_indices != metadata['mtf_result']:
                print(f"WARNING: Decoded MTF indices don't match stored values")
        except ValueError as e:
            raise ValueError(f"Huffman decoding failed: {str(e)}")
        
        # Apply inverse MTF
        print("Decoding Move-to-Front...")
        alphabet = metadata['alphabet']
        bwt_result = MoveToFront.decode(mtf_indices, alphabet)
        print(f"Decoded BWT result: '{bwt_result[:20]}{'...' if len(bwt_result) > 20 else ''}'")
        
        # Verify BWT result if available for debugging
        if 'bwt_result' in metadata and bwt_result != metadata['bwt_result']:
            print(f"WARNING: Decoded BWT result doesn't match expected value")
        
        # Apply inverse BWT
        print("Applying inverse BWT...")
        original_seq = BurrowsWheeler.inverse_bwt(bwt_result)
        print(f"Recovered sequence: '{original_seq[:20]}{'...' if len(original_seq) > 20 else ''}'")
        
        # Verify length matches original
        if len(original_seq) != metadata['original_length']:
            print(f"Warning: Length mismatch! Original: {metadata['original_length']}, Decompressed: {len(original_seq)}")
        
        # Write output
        print(f"Writing output to {output_path}...")
        with open(output_path, 'w') as file:
            file.write(original_seq)
        
        print(f"Decompression complete. Output saved to {output_path}")
        print(f"Decompressed size: {len(original_seq)} bytes")
        return output_path

def main():
    parser = argparse.ArgumentParser(
        description='DNA Sequence Compressor using BWT, Move-to-Front, and Huffman Coding',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Compression command
    compress_parser = subparsers.add_parser('compress', help='Compress a DNA sequence file')
    compress_parser.add_argument('input', help='Input FASTA file')
    compress_parser.add_argument('-o', '--output', help='Output compressed file (default: input.compressed)')
    compress_parser.add_argument('-c', '--chunk', type=int, default=1000000, 
                               help='Chunk size in bytes (for future use)')
    compress_parser.add_argument('--no-suffix-array', action='store_true',
                               help='Disable suffix array optimization for BWT (use original algorithm)')

    # Decompression command
    decompress_parser = subparsers.add_parser('decompress', help='Decompress a compressed file')
    decompress_parser.add_argument('input', help='Input compressed file (without .meta extension)')
    decompress_parser.add_argument('-o', '--output', help='Output FASTA file (default: input.decompressed.fasta)')

    args = parser.parse_args()

    try:
        if args.command == 'compress':
            output = args.output if args.output else args.input + '.compressed'
            compressor = DNACompressor(args.input, args.chunk, not args.no_suffix_array)
            compressor.compress(output)
        elif args.command == 'decompress':
            output = args.output if args.output else args.input.replace('.compressed', '') + '.decompressed.fasta'
            compressor = DNACompressor(args.input)
            compressor.decompress(args.input, output)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()