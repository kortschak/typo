// Copyright ©2017 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package typo implements the game of typogenetics as described in Gödel, Escher
// Bach: an Eternal Golden Braid by Douglas Hofstadter.
package typo

import (
	"bytes"
	"fmt"
	"text/tabwriter"
)

// AminoAcid represents a typogenetics amino acid.
type AminoAcid byte

const (
	Non AminoAcid = iota

	Cut // cut strand(s)
	Del // delete a base from strand
	Swi // switch enzyme to other strand
	Mvr // move one unit to the right
	Mvl // move one unit to the left
	Cop // turn on Copy mode
	Off // turn off Copy mode
	Ina // insert A to the right of this unit
	Inc // insert C to the right of this unit
	Ing // insert G to the right of this unit
	Int // insert T to the right of this unit
	Rpy // search for the nearest pyrimidine to the right
	Rpu // search for the nearest purine to the right
	Lpy // search for the nearest pyrimidine to the left
	Lpu // search for the nearest purine to the left
)

func (a AminoAcid) String() string {
	return [...]string{
		"Non", "Cut", "Del", "Swi",
		"Mvr", "Mvl", "Cop", "Off",
		"Ina", "Inc", "Ing", "Int",
		"Rpy", "Rpu", "Lpy", "Lpu",
	}[a]
}

// Direction represent the tertiary structure direction of a typogenetics enzyme.
type Direction int8

const (
	North Direction = iota
	East
	South
	West
)

func (d Direction) String() string {
	return [...]string{"north", "east", "south", "west"}[d]
}

// Kink represents an enzyme folding operation.
type Kink int8

const (
	Left = iota - 1
	Straight
	Right
)

var (
	// Complement describes how bases pair.
	Complement = [255]func() byte{'A': T, 'C': G, 'G': C, 'T': A}

	// Code is the pseudo-codon to amino acid lookup.
	Code = [16]AminoAcid{
		Non, Cut, Del, Swi,
		Mvr, Mvl, Cop, Off,
		Ina, Inc, Ing, Int,
		Rpy, Rpu, Lpy, Lpu,
	}
	// index is a lookup into the Code table.
	index = [255]byte{'A': 0, 'C': 1, 'G': 2, 'T': 3}

	// Inserts specifies what base or pseudobase
	// is inserted by insert and cut amino acids.
	Inserts = [16]func() byte{
		Cut: null, // This value must not be altered.
		Ina: A,
		Inc: C,
		Ing: G,
		Int: T,
	}

	// Moves specifies the direction the search amino
	// acids move.
	Moves = [16]int{
		Rpy: Right,
		Rpu: Right,
		Lpy: Left,
		Lpu: Left,
	}

	// Matches specifies the matching criteria of
	// search amino acids.
	Matches = [16]func(byte) bool{
		Rpy: IsPyrimidine,
		Rpu: IsPurine,
		Lpy: IsPyrimidine,
		Lpu: IsPurine,
	}

	// Kinks specifies the folding characteristics
	// of each amino acid.
	Kinks = [16]Kink{
		Cut: Straight,
		Del: Straight,
		Swi: Right,
		Mvr: Straight,
		Mvl: Straight,
		Cop: Right,
		Off: Left,
		Ina: Straight,
		Inc: Right,
		Ing: Right,
		Int: Left,
		Rpy: Right,
		Rpu: Left,
		Lpy: Left,
		Lpu: Left,
	}

	// Preference specifies the binding preference of
	// each folding direction.
	Preference = [4]func() byte{
		East:  A,
		North: C,
		South: G,
		West:  T,
	}
)

func init() {
	for b := range index {
		switch b {
		default:
			index[b] = 0xff
		case 'A', 'C', 'G', 'T':
		}
	}
}

func null() byte { return 0 }

// A returns the 'A' base.
func A() byte { return 'A' }

// C returns the 'C' base.
func C() byte { return 'C' }

// G returns the 'G' base.
func G() byte { return 'G' }

// T returns the 'T' base.
func T() byte { return 'T' }

// IsPurine returns whether b is 'A' or 'G'.
func IsPurine(b byte) bool { return b == 'A' || b == 'G' }

// IsPurine returns whether b is 'C' or 'T'.
func IsPyrimidine(b byte) bool { return b == 'C' || b == 'T' }

// Enzyme is implements a typogenetics enzyme.
type Enzyme []AminoAcid

func (e Enzyme) String() string {
	var buf bytes.Buffer
	for i, c := range e {
		if i != 0 {
			fmt.Fprint(&buf, "-")
		}
		fmt.Fprint(&buf, c)
	}
	return buf.String()
}

// Fold returns the first and last segment folding directions of the receiver.
func (e Enzyme) Fold() (first, last Direction) {
	dir := North
	// The text is ambiguous about the behaviour here; the folding
	// example on p511 gives Rpy-Ina-Rpu-Mvr-Int-Mvl-Cut-Swi-Cop
	// as folding with ⇒ ⇑, though Cop is a right-turning amino
	// acid. The text states that the segments are perpendicular,
	// but if we count all segments, they should be parallel (⇒⇒).
	// The only other example is the operation example on p508
	// with the enzyme Rpu-Inc-Cop-Mvr-Mvl-Swi-Lpu-Int which Hofstadter
	// says has a preference for 'G', meaning it must have a ⇒⇓
	// fold. The only way to obtain this result is to count all
	// segments.
	for _, d := range e {
		dir += Direction(Kinks[d])
		dir &= 0x3
	}
	return Direction(Kinks[e[0]] & 0x3), dir
}

// Preference returns the base preference of the receiver.
func (e Enzyme) Preference() byte {
	first, last := e.Fold()
	t := Kink(East - first)
	return Preference[(last+Direction(t))&0x3]()
}

// OperateOn performs the enzymatic activity of the receiver on the given
// typogenetics complex starting from the specified position of the first
// strand of the complex according to rules of typogenetics on pp504-513
// of GEB and returns the resulting product complex.
// If debug is not nil, the sequence of operations and the intermediate
// results are written into the buffer.
// OperateOn will panic if the receiver includes an unknown amino acid
// or the Non amino acid.
func (e Enzyme) OperateOn(c Complex, pos int, debug *bytes.Buffer) Complex {
	copyMode := false
	if len(c[0]) != len(c[1]) {
		panic("typo: invalid Complex: length mismatch")
	}
	var w *tabwriter.Writer
	if debug != nil {
		w = tabwriter.NewWriter(debug, 0, 0, 1, ' ', 0)
	}
	completed := true
	for _, cmd := range e {
		if w != nil {
			fmt.Fprintf(w, "%s\t%4d\t%c\t%q\t%q\n", cmd, pos, c[0][pos], c[0], c[1])
		}

		switch cmd {
		default:
			panic("unknown amino acid in enzyme")
		case Non:
			panic("non used in enzyme")
		case Del:
			c[0][pos] = 0
			pos--
		case Swi:
			c[0], c[1] = c[1], c[0]
			pos = len(c[0]) - pos - 1
		case Mvr:
			pos++
		case Mvl:
			pos--
		case Cop:
			copyMode = true
		case Off:
			copyMode = false
		case Cut, Ina, Inc, Ing, Int:
			pos++
			c[0] = insert(c[0], Inserts[cmd](), pos, false)
			c[1] = insert(c[1], 0, pos, true)
		case Rpy, Rpu, Lpy, Lpu:
			move := Moves[cmd]
			pos += move
			for found := Matches[cmd]; onStrand(c[0], pos); pos += move {
				if copyMode {
					copyOpposite(c[1], c[0], pos)
				}
				if found(c[0][pos]) {
					break
				}
			}
		}
		if !onStrand(c[0], pos) {
			if w != nil {
				if 0 <= pos && pos < len(c[0]) {
					fmt.Fprintf(w, "empty\t%4d\t·\t%q\t%q\n", pos, c[0], c[1])
				} else {
					fmt.Fprintf(w, "off\t%4d\t-\t%q\t%q\n", pos, c[0], c[1])
				}
			}
			completed = false
			break
		}
		if copyMode {
			copyOpposite(c[1], c[0], pos)
		}
	}
	if w != nil {
		if completed {
			var b byte
			if 0 <= pos && pos < len(c[0]) {
				b = c[0][pos]
				if b == 0 {
					b = '·'
				}
			} else {
				b = '-'
			}
			fmt.Fprintf(w, "done\t%4d\t%c\t%q\t%q\n", pos, b, c[0], c[1])
		}
		w.Flush()
	}

	return c
}

func onStrand(s Strand, pos int) bool {
	return 0 <= pos && pos < len(s) && s[pos] != 0
}

// Strand is a typogenetics base strand.
type Strand []byte

func (s Strand) String() string {
	var buf bytes.Buffer
	for _, b := range s {
		if b == 0 {
			buf.WriteRune('·')
			continue
		}
		buf.WriteByte(b)
	}
	return buf.String()
}

func copyOpposite(dst, src Strand, pos int) {
	if index[src[pos]] != 0xff {
		dst[len(dst)-pos-1] = Complement[src[pos]]()
	}
}

func prepend(s Strand, b ...byte) Strand {
	return append(b, s...)
}

func insert(s Strand, b byte, pos int, opposite bool) Strand {
	if opposite {
		pos = len(s) - pos
	}
	if pos == len(s) {
		return append(s, b)
	}
	if pos == -1 {
		return prepend(s, b)
	}
	s = append(s[:pos+1], s[pos:]...)
	s[pos] = b
	return s
}

// Enzymes returns the set of enzymes specified by the receiver.
func (s Strand) Enzymes() []Enzyme {
	buf := make([]AminoAcid, 0, len(s)/2)
	var e []Enzyme
	for i := 0; i+1 < len(s); i += 2 {
		c := Code[index[s[i]]*4+index[s[i+1]]]
		if c == Non {
			if len(buf) != 0 {
				e = append(e, buf)
				buf = buf[len(buf):]
			}
			continue
		}
		buf = append(buf, c)
	}
	if len(buf) != 0 {
		e = append(e, buf)
	}
	return e
}

// Complex is a complex of complementary strands.
type Complex [2]Strand

// NewComplex returns a new valid complex.
func NewComplex(s Strand) Complex { return Complex{s, make(Strand, len(s))} }

// Products returns the dissociated strands of a complex.
func (c Complex) Products() []Strand {
	var n int
	for _, s := range c {
		var isSeq bool
		for _, b := range s {
			wasSeq := isSeq
			isSeq = b != 0
			if isSeq && !wasSeq {
				n++
			}
		}
	}

	p := make([]Strand, 0, n)
	for _, s := range c {
		var isSeq bool
		start := -1
		for i, b := range s {
			wasSeq := isSeq
			isSeq = b != 0
			if isSeq && !wasSeq {
				start = i
			} else if wasSeq && !isSeq {
				p = append(p, s[start:i])
				start = -1
			}
		}
		if isSeq {
			p = append(p, s[start:])
		}
	}

	return p
}
