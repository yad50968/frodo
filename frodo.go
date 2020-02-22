package frodo

import "encoding/binary"

// Frodo interface
type Frodo interface {
	Encode(k []byte) [][]uint16                   // Encode encodes an integer 0 ≤ k < 2^B as an element in Zq by multiplying it by q/2B = 2^(D−B): ec(k) := k·q/2^B
	Decode(K [][]uint16) []byte                   // Decode decodes the m-by-n matrix K into a bit string of length l = B·m·n. dc(c) = ⌊c·2^B/q⌉ mod 2^B
	Pack(C [][]uint16) []byte                     // Pack packs a matrix into a bit string
	Unpack(b []byte, n1, n2 int) [][]uint16       // Unpack unpacks a bit string into a matrix n1-by-n2
	Gen(seed []byte) [][]uint16                   // Gen returns a pseudorandom matrix using SHAKE128
	Sample(t uint16) uint16                       // Sample returns a sample e from the distribution χ
	SampleMatrix(r []byte, n1, n2 int) [][]uint16 // SampleMatrix sample the n1-by-n2 matrix entry
}

// Parameters of frodo KEM mechanism
type Parameters struct {
	no     int    // n ≡ 0 (mod 8) the main parameter
	q      uint16 // a power-of-two integer modulus with exponent D ≤ 16 !! minus one for bit masking
	D      int    // a power
	m, n   int    // integer matrix dimensions with
	B      int    // the number of bits encoded in each matrix entry
	l      int    // B·m·n, the length of bit strings that are encoded as m-by-n matrices
	lseedA int    // the byte length of public matrix A
}

// Frodo640 returns Parameters struct no.640
func Frodo640() *Parameters {

	param := new(Parameters)

	param.no = 640
	param.q = 0x7fff
	param.D = 15
	param.B = 2
	param.m = 8
	param.n = 8
	param.lseedA = 16
	param.l = 16

	return param
}

// Frodo976 returns Parameters struct no.976
func Frodo976() *Parameters {

	param := new(Parameters)

	param.no = 976
	param.q = 0xffff
	param.D = 16
	param.B = 3
	param.m = 8
	param.n = 8
	param.lseedA = 16
	param.l = 24

	return param
}

// Frodo1344 returns Parameters struct no.1344
func Frodo1344() *Parameters {

	param := new(Parameters)

	param.no = 1344
	param.q = 0xffff
	param.D = 16
	param.B = 4
	param.m = 8
	param.n = 8
	param.lseedA = 16
	param.l = 32

	return param
}

// Encode encodes an integer 0 ≤ k < 2^B as an element in Zq
// by multiplying it by q/2B = 2^(D−B): ec(k) := k·q/2^B
func (param *Parameters) Encode(k []byte) [][]uint16 {

	K := make([][]uint16, param.m)
	for i := range K {
		K[i] = make([]uint16, param.n)
		for j := range K[i] {
			temp := uint16(0)
			for l := 0; l < param.B; l++ {
				index, shift := ((i*param.n+j)*param.B+l)/8, uint(((i*param.n+j)*param.B+l)&7)
				if k[index]&(byte(0x80)>>shift) != 0 { // litte-endian
					temp |= uint16(1 << uint(l))
				}
			}
			K[i][j] = param.ec(temp)
		}
	}
	return K
}

// Decode decodes the m-by-n matrix K into a bit string of {0,1}^(B·m·n). dc(c) = ⌊c·2^B/q⌉ mod 2^B
func (param *Parameters) Decode(K [][]uint16) []byte {

	k := make([]byte, param.l)
	for i, row := range K {
		for j := range row {
			temp := param.dc(K[i][j])
			for l := 0; l < param.B; l++ {
				if temp&uint16(1<<uint(l)) != 0 {
					index, shift := ((i*param.n+j)*param.B+l)/8, uint(((i*param.n+j)*param.B+l)&7)
					k[index] |= byte(0x80) >> shift // litte-endian
				}
			}
		}
	}
	return k
}

// Pack packs a n1-by-n2 matrix over Zq into a bit string {0,1}^(D*n1*n2)
func (param *Parameters) Pack(C [][]uint16) []byte {

	n1, n2 := len(C), len(C[0])
	b := make([]byte, param.D*n1*n2/8)
	for i := 0; i < n1; i++ {
		for j := 0; j < n2; j++ {
			for l := 0; l < param.D; l++ {
				if (uint16(1)<<uint(param.D-1-l))&C[i][j] != 0 {
					index, shift := ((i*n2+j)*param.D+l)/8, uint(((i*n2+j)*param.D+l)&7)
					b[index] |= byte(0x80) >> shift
				}
			}
		}
	}
	return b
}

// Unpack unpacks a bit string {0,1}^(D*n1*n2) into a matrix (n1-by-n2) over Zq
func (param *Parameters) Unpack(b []byte, n1, n2 int) [][]uint16 {

	C := make([][]uint16, n1)
	for i := range C {
		C[i] = make([]uint16, n2)
		for j := range C[i] {
			for l := 0; l < param.D; l++ {
				index, shift := ((i*n2+j)*param.D+l)/8, uint(((i*n2+j)*param.D+l)&7)
				if b[index]&byte(0x80>>shift) != 0 {
					C[i][j] |= uint16(1) << uint(param.D-1-l)
				}
			}
		}
	}
	return C
}

// Gen returns a pseudorandom matrix using SHAKE128
func (param *Parameters) Gen(seed []byte) [][]uint16 {

	A := make([][]uint16, param.no)
	for i := uint16(0); i < uint16(param.no); i++ {

		b := []byte{byte(i >> 8), byte(i)}
		b = append(b, seed...)
		shakeStr := param.shake(b, param.no*2)

		A[i] = make([]uint16, param.no)
		for j := 0; j < param.no; j++ {
			A[i][j] = ((uint16(shakeStr[j*2]) << 8) | uint16(shakeStr[i*2+1])) & param.q
		}
	}

	return A
}

// SampleMatrix sample the n1-by-n2 matrix entry
func (param *Parameters) SampleMatrix(r []byte, n1, n2 int) [][]uint16 {

	E := make([][]uint16, n1)
	for i := 0; i < n1; i++ {
		E[i] = make([]uint16, n2)
		for j := 0; j < n2; j++ {
			t, d := binary.LittleEndian.Uint16(r[2*(i*n2+j):]), uint32(0)
			r := (uint32(t) << 8) + uint32(t)
			for j := uint(0); j < 8; j++ {
				d += (r << j) & 0x01010101
			}
			a := ((d >> 8) & 0xff) + (d & 0xff)
			b := (d >> 24) + ((d >> 16) & 0xff)
			E[i][j] = 0xfffd - uint16(a) + uint16(b)
		}
	}

	return E
}
