for ((st=4401;st<64000;st+=20))
   do
   python /media/KINGSTON/current-2-2-11/maple/zeros.py foo.txt $st &
   st="`expr $st + 20`"
   python /media/KINGSTON/current-2-2-11/maple/zeros.py bar.txt $st &
   st="`expr $st + 20`"
   python /media/KINGSTON/current-2-2-11/maple/zeros.py bletch.txt $st &
   wait
   cat foo.txt bar.txt bletch.txt >> zeros0.txt
done