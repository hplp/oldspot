TARGET=oldspot

SRCDIR=src
SRC=$(wildcard $(SRCDIR)/*.cc)
INCDIR=src
OBJDIR=obj
OBJ=$(SRC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

OPT=-O3
INCLUDE=-I$(INCDIR)
CXXFLAGS += -std=c++11 -Wall $(INCLUDE) $(OPT)
LIBS=-lm -lpugixml
LFLAGS += $(LIBS) $(OPT)

.PHONY: $(TARGET) debug clean

default: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $^ -o $@ $(LFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@mkdir -p $(OBJDIR)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

debug: CXXFLAGS += -g
debug: LFLAGS += -g
debug: OPT=-O0
debug: $(TARGET)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(TARGET)