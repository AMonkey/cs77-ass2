#!/bin/bash

function check_error {
    if [ "$?" -ne "0" ]; then
        echo "An error occurred!"
        exit 1
    fi
}

rm -f scenes/*.png

echo "test01: Cylinder"
./view-ref-linux -i scenes/test01.json
check_error
mv scenes/test01.png scenes/test01-ref.png

echo "test02: Mesh with triangles"
./view-ref-linux -i scenes/test02.json
check_error
mv scenes/test02.png scenes/test02-ref.png

echo "test03: Mesh with quads and triangles"
./view-ref-linux -i scenes/test03.json
check_error
mv scenes/test03.png scenes/test03-ref.png

echo "test04: Catmull-Clark Subdiv"
./view-ref-linux -i scenes/test04.json
check_error
mv scenes/test04.png scenes/test04-ref.png

echo "test05: Transforms"
./view-ref-linux -i scenes/test05.json
check_error
mv scenes/test05.png scenes/test05-ref.png

echo "test06: Bezier Spline"
./view-ref-linux -i scenes/test06.json
check_error
mv scenes/test06.png scenes/test06-ref.png

./render_all.sh
