if [ -z "$1" ]; then
  echo "Usage: $0 <filename>"
  exit 1
fi

filename="${1%.*}"

g++ "$filename.cpp" -fopenmp -o "$filename"

if [ $? -eq 0 ]; then
  ./"$filename"
else
  echo "Compilation failed."
  exit 1
fi