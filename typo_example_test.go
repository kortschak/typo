// Copyright ©2017 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package typo

import (
	"bytes"
	"fmt"
)

func Example_book() {
	e := Enzyme{Rpu, Inc, Cop, Mvr, Mvl, Swi, Lpu, Int}
	s := "TAGATCCAGTCCATCGA"

	first, last := e.Fold()
	pref := e.Preference()
	fmt.Printf("Folding:\n First:%s Last:%s\n\nPrefers:%q\n\n", first, last, pref)

	var pos []int
	for i, b := range Strand(s) {
		if b == pref {
			pos = append(pos, i)
		}
	}

	fmt.Printf("%d possible start positions: %d\n\n", len(pos), pos)
	for _, p := range pos {
		var buf bytes.Buffer
		fmt.Printf("Start at %d:\n\n", p)
		// As stated in the game rules, the strand
		// being operated on is consumed by the enzyme
		// so this example creates a new Strand from
		// a string constate for each start position.
		products := e.OperateOn(NewComplex(Strand(s)), p, &buf).Products()
		fmt.Printf("%s\nProducts:%q\n\n", &buf, products)
	}

	// Output:
	//
	// Folding:
	//  First:west Last:north
	//
	// Prefers:'G'
	//
	// 3 possible start positions: [2 8 15]
	//
	// Start at 2:
	//
	// Rpu     2 G "TAGATCCAGTCCATCGA"   "·················"
	// Inc     3 A "TAGATCCAGTCCATCGA"   "·················"
	// Cop     4 C "TAGACTCCAGTCCATCGA"  "··················"
	// Mvr     4 C "TAGACTCCAGTCCATCGA"  "·············G····"
	// Mvl     5 T "TAGACTCCAGTCCATCGA"  "············AG····"
	// Swi     4 C "TAGACTCCAGTCCATCGA"  "············AG····"
	// Lpu    13 G "············AG····"  "TAGACTCCAGTCCATCGA"
	// Int    12 A "············AG····"  "TAGACTCCAGTCCATCGA"
	// done   13 T "············ATG····" "TAGACATCCAGTCCATCGA"
	//
	// Products:["ATG" "TAGACATCCAGTCCATCGA"]
	//
	// Start at 8:
	//
	// Rpu     8 G "TAGATCCAGTCCATCGA"   "·················"
	// Inc    12 A "TAGATCCAGTCCATCGA"   "·················"
	// Cop    13 C "TAGATCCAGTCCACTCGA"  "··················"
	// Mvr    13 C "TAGATCCAGTCCACTCGA"  "····G·············"
	// Mvl    14 T "TAGATCCAGTCCACTCGA"  "···AG·············"
	// Swi    13 C "TAGATCCAGTCCACTCGA"  "···AG·············"
	// Lpu     4 G "···AG·············"  "TAGATCCAGTCCACTCGA"
	// Int     3 A "···AG·············"  "TAGATCCAGTCCACTCGA"
	// done    4 T "···ATG·············" "TAGATCCAGTCCACATCGA"
	//
	// Products:["ATG" "TAGATCCAGTCCACATCGA"]
	//
	// Start at 15:
	//
	// Rpu   15 G "TAGATCCAGTCCATCGA"  "·················"
	// Inc   16 A "TAGATCCAGTCCATCGA"  "·················"
	// Cop   17 C "TAGATCCAGTCCATCGAC" "··················"
	// Mvr   17 C "TAGATCCAGTCCATCGAC" "G·················"
	// off   18 - "TAGATCCAGTCCATCGAC" "G·················"
	//
	// Products:["TAGATCCAGTCCATCGAC" "G"]
}

func Example_quine() {
	seed := "CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG"
	pool := []Strand{Strand(seed)}

	fmt.Printf("%s: 1\n", seed)

	for n := 0; n < 5; n++ {
		var daughters []Strand
		for _, s := range pool {
			m := s.Enzymes()
			c := NewComplex(s)
			for _, e := range m {
				pref := e.Preference()
				for i, b := range c[0] {
					if b == pref {
						c = e.OperateOn(c, i, nil)
						break
					}
				}
			}
			daughters = append(daughters, c.Products()...)
		}
		pool = daughters

		gen := make(map[string]int)
		for _, s := range pool {
			gen[string(s)]++
		}
		for s, count := range gen {
			fmt.Printf("%s: %d\n", s, count)
		}
	}

	// Output:
	//
	// CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG: 1
	// CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG: 2
	// CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG: 4
	// CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG: 8
	// CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG: 16
	// CGTTCCTCTCTCTCTATAGAGAGAGAGGAACG: 32
}
