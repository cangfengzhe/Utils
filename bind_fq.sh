#!/usr/bin/env bash
top_dir=/leostore/Demultiplex/xten/x10samba495/170123_E00495_0091_AHHJN5ALXX/Demultiplex_L8_8/HHJN5ALXX
bottom_dir=/leostore/Demultiplex/xten/x10samba495/170123_E00495_0091_AHHJN5ALXX/Demultiplex_L7_8/HHJN5ALXX
output_dir=/leostore/Demultiplex/xten/x10samba495/170123_E00495_0091_AHHJN5ALXX/Run722_bind
n=0
for xx in `ls $top_dir`;do
    cat $top_dir/$xx $bottom_dir/$xx > $output_dir/$xx &
    n=$(($n+1))    
    q=$(( $n % 5 ))
    #echo $q
    if [ $q = 0 ];then
        wait
    fi

    
done
