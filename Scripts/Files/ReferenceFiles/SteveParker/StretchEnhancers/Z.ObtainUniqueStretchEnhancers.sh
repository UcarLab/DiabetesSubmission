#!/bin/bash

cat /Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/TissuesStretchEnhancers.txt | while read LINE
do
	cat $LINE > temp.bed
	cat /Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/TissuesStretchEnhancers.txt | while read SENTENCE
	do
		if [[ $LINE = $SENTENCE ]]
		then
			continue
		else
			intersectBed -a temp.bed -b $SENTENCE -v | uniq > temp2.bed
			mv temp2.bed temp.bed
		fi
	done
	cat temp.bed | sort -k1,1 -k2,2n | uniq > "/Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/StretchEnhancers/Z.Unique_StretchEnhancers/"$LINE
done

cat /Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/StretchEnhancers/Z.Unique_StretchEnhancers/*.bed | sort -k1,1 -k2,2n | uniq | bedtools merge -i - | sort -k1,1 -k2,2n | uniq > /Users/skn/Desktop/Lab/Diabetic_Versus_Non-Diabetic/ReferenceFiles/SteveParker/StretchEnhancers/Z.Unique_StretchEnhancers/0.All_StretchEnhancers_Merged.bed
