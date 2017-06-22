// Copyright ©2017 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package typo

import (
	"bytes"
	"reflect"
	"testing"
)

var foldingTests = []struct {
	enzyme Enzyme
	want   Direction
}{
	{enzyme: Enzyme{Rpy, Ina, Rpu, Mvr, Int, Mvl, Cut, Swi, Cop}, want: East},
	{enzyme: Enzyme{Rpu, Inc, Cop, Mvr, Mvl, Swi, Lpu, Int}, want: North},
	{enzyme: Enzyme{Ina, Rpu, Cop, Inc, Swi}, want: South},
}

func TestEnzymeFold(t *testing.T) {
	for _, test := range foldingTests {
		_, got := test.enzyme.Fold()
		if got != test.want {
			t.Errorf("unexpected fold direction: got:%q want:%q", got, test.want)
		}
	}
}

var preferenceTests = []struct {
	enzyme Enzyme
	want   byte
}{
	{enzyme: Enzyme{Rpy, Ina, Rpu, Mvr, Int, Mvl, Cut, Swi, Cop}, want: 'A'},
	{enzyme: Enzyme{Rpu, Inc, Cop, Mvr, Mvl, Swi, Lpu, Int}, want: 'G'},
	{enzyme: Enzyme{Ina, Rpu, Cop, Inc, Swi}, want: 'T'},
}

func TestEnzymePreference(t *testing.T) {
	for _, test := range preferenceTests {
		got := test.enzyme.Preference()
		if got != test.want {
			t.Errorf("unexpected preference: got:%q want:%q", got, test.want)
		}
	}
}

var strandEnzymeTests = []struct {
	strand Strand
	want   []Enzyme
}{
	{strand: []byte("TAGATCCAGTCCACATCGA"), want: []Enzyme{{Rpy, Ina, Rpu, Mvr, Int, Mvl, Cut, Swi, Cop}}},
	{strand: []byte("GATCCGGCAT"), want: []Enzyme{{Ina, Rpu, Cop, Inc, Swi}}},
	{strand: []byte("AA"), want: nil},
	{strand: []byte("AAAA"), want: nil},
}

func TestStrandEnzyme(t *testing.T) {
	for _, test := range strandEnzymeTests {
		got := test.strand.Enzymes()
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("unexpected preference: got:%q want:%q", got, test.want)
		}
	}
}

var copyOppositeTests = []struct {
	strand Strand
	pos    int
	want   Strand
}{
	{strand: []byte("xAGATCCAGTCCACATCGA"), pos: 0, want: []byte("                   ")},
	{strand: []byte("TAGATCCAGTCCACATCGA"), pos: 0, want: []byte("                  A")},
	{strand: []byte("TAGATCCAGTCCACATCGA"), pos: 1, want: []byte("                 T ")},
	{strand: []byte("TAGATCCAGTCCACATCGA"), pos: 2, want: []byte("                C  ")},
	{strand: []byte("TAGATCCAGTCCACATCGA"), pos: 16, want: []byte("  G                ")},
	{strand: []byte("TAGATCCAGTCCACATCGA"), pos: 17, want: []byte(" C                 ")},
	{strand: []byte("TAGATCCAGTCCACATCGA"), pos: 18, want: []byte("T                  ")},
}

func TestCopyOpposite(t *testing.T) {
	for _, test := range copyOppositeTests {
		got := make(Strand, len(test.strand))
		for i := range got {
			got[i] = ' '
		}
		copyOpposite(got, test.strand, test.pos)
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("unexpected prepend result:\ngot: %q\nwant:%q", got, test.want)
		}
	}
}

var prependTests = []struct {
	strand Strand
	add    byte
	want   Strand
}{
	{strand: []byte("TAGATCCAGTCCACATCGA"), add: 'C', want: []byte("CTAGATCCAGTCCACATCGA")},
}

func TestPrepend(t *testing.T) {
	for _, test := range prependTests {
		got := prepend(test.strand, test.add)
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("unexpected prepend result:\ngot: %q\nwant:%q", got, test.want)
		}
	}
}

var insertTests = []struct {
	strand   Strand
	add      byte
	pos      int
	opposite bool
	want     Strand
}{
	{strand: []byte("A"), add: 'G', pos: 1, opposite: false, want: []byte("AG")},
	{strand: []byte("A"), add: 'G', pos: 1, opposite: true, want: []byte("GA")},
	{strand: []byte("ACTGC"), add: 'G', pos: 2, opposite: false, want: []byte("ACGTGC")},
	{strand: []byte("GCAGT"), add: ' ', pos: 2, opposite: true, want: []byte("GCA GT")},
	{strand: []byte("TAGATCCAGTCCATCGA"), add: 'C', pos: 13, want: []byte("TAGATCCAGTCCACTCGA")},
}

func TestInsert(t *testing.T) {
	for _, test := range insertTests {
		got := insert(test.strand, test.add, test.pos, test.opposite)
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("unexpected insert result:\ngot: %q\nwant:%q", got, test.want)
		}
	}
}

var productsTests = []struct {
	c    Complex
	want []Strand
}{
	{c: Complex{{}, {}}, want: []Strand{}},
	{
		c: Complex{
			{0, 0, 1, 2, 3, 0, 0, 4, 0, 5},
			{0, 6, 7, 0, 8, 0, 0, 9, 10, 0},
		},
		want: []Strand{
			{1, 2, 3},
			{4},
			{5},
			{6, 7},
			{8},
			{9, 10},
		},
	},
	{
		c: Complex{
			{1, 2, 3, 0, 0, 4, 0, 5},
			{6, 7, 8, 0, 0, 9, 10, 0},
		},
		want: []Strand{
			{1, 2, 3},
			{4},
			{5},
			{6, 7, 8},
			{9, 10},
		},
	},
}

func TestProducts(t *testing.T) {
	for _, test := range productsTests {
		got := test.c.Products()
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("unexpected products result:\ngot: %q\nwant:%q", got, test.want)
		}
	}
}

var operateOnTests = []struct {
	enzyme Enzyme
	strand Strand
	pos    int
	want   []Strand
}{
	{
		enzyme: Enzyme{Rpu, Inc, Cop, Mvr, Mvl, Swi, Lpu, Int},
		strand: Strand("TAGATCCAGTCCATCGA"),
		pos:    8,
		want: []Strand{
			Strand("ATG"),
			Strand("TAGATCCAGTCCACATCGA"),
		},
	},
	{
		enzyme: Enzyme{Rpu, Inc, Cop, Mvr, Mvl, Swi, Lpu, Int},
		strand: Strand("TAGATCCAGTCCATCGA"),
		pos:    15,
		want: []Strand{
			Strand("TAGATCCAGTCCATCGAC"),
			Strand("G"),
		},
	},
	{
		enzyme: Enzyme{Ina, Rpu, Cop, Inc, Swi},
		strand: Strand("GATCCGGCAT"),
		pos:    2,
		want: []Strand{
			Strand("GC"),
			Strand("GATACCGCGCAT"),
		},
	},
}

func TestOperateOn(t *testing.T) {
	for _, test := range operateOnTests {
		var buf bytes.Buffer
		got := test.enzyme.OperateOn(NewComplex(test.strand), test.pos, &buf).Products()
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("unexpected operation results:\ngot: %q\nwant:%q", got, test.want)
			t.Logf("\n%s", &buf)
		}
	}
}
