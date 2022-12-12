while true; do
    file=$(inotifywait /home/Mathieu/Documents/TRAPPIST/tmpout/ -q -e modify,create --format %f)

    if [[ $file == *".png"* ]];
    then
        kitty +kitten icat $file &
    fi
done
