import sys

import filter_mtz


def main():
  filter_mtz.run(sys.argv[1:])

if __name__ == '__main__':
  main()
