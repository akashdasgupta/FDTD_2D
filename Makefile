CC = mpic++
IDIR = include

CFLAGS = -I $(IDIR) -O2

SRC = $(wildcard *.cpp)
NAME = o2_sweep_exec2


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC) -Wno-narrowing

.PHONY: clean
clean:
	@rm -f $(NAME) data/* runfile.* movie_runfile.* *.mp4 
