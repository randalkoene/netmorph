#! /bin/sh

#imgnum=9999
imgnum=1000

for n in `ls $1*.fig`; do

    echo $n
    fig2dev -L png -Z$2 $n $1-$imgnum.png

    #imgnum=`expr $imgnum - 1`
    imgnum=`expr $imgnum + 1`

done

ls $1-*.png
printf "Remove the .fig sample files (y/N)?"
read rmfigs

if [ "$rmfigs" = "y" ]; then
    rm -f $1*.fig
fi

#imgnum=`expr $imgnum + 1`
imgnum=1000

gimp $1-$imgnum.png

ls $1*gif
printf "Remove the .png sample files (y/N)?"
read rmfigs

if [ "$rmfigs" = "y" ]; then
    rm -f $1*.png
fi
