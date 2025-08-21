#!/bin/bash

for i in {1..12}
do
	mkdir $i
	cp b.in $i/b.in
	cp c.in $i/c.in
	sed -i "s/{lambdaX}/{lambda$i}/g" $i/b.in
	sed -i "s/{lambdaX}/{lambda$i}/g" $i/c.in
	cp lmp.data $i/lmp.data
	cp lmp.setting $i/lmp.setting
done