import sys

from . import slicemtz


def main():
  slicemtz.run(sys.argv[1:])

if __name__ == '__main__':
  main()
