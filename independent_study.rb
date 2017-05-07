#!/usr/bin/env jruby
require 'benchmark'
require 'set'
require 'java'
require '/Users/YASH/RubymineProjects/IndependentStudy/Matlab4.jar'

class Param

  MIN_INTERVAL = 60.0
  MAX_INTERVAL = 100.0
  MIN_PROB = 0.5
  MAX_PROB = 0.6
  MIN_COST = 0.5
  MAX_COST = 0.6

end
# noinspection ALL
class Assign
  @@assigner = nil;
  @@remainingHelpers = [];
  @@helperList = [];

  #same initial solutions for simulated annealing, best local search and multiple start BLS
  @@initialSolutionBLS = nil
  @@initialCostBLS = nil

  #create an assigner object
  def Assign.create
    @@assigner = new unless @@assigner
    @@assigner
  end

  #expression = @@clause, helperList = @@helpers
  def init(numberOfSensors)
    tau = 0.5
    duration = 10000
    lmd = Array.new(numberOfSensors)
    p = Array.new(numberOfSensors)
    c = Array.new(numberOfSensors)
    c_ntfy = Array.new(numberOfSensors)
    count = 0
    numberOfSensors.times do
      p0 = Param::MIN_PROB + Random.rand(0..(Param::MAX_PROB - Param::MIN_PROB))
      p[count] = p0
      intv = Param::MIN_INTERVAL + Random.rand(0..(Param::MAX_INTERVAL - Param::MIN_INTERVAL))
      ld = Float(1)/intv
      lmd[count] = ld
      cread = Param::MIN_COST + Random.rand(0..(Param::MAX_COST - Param::MIN_COST))
      c[count] = cread
      c_ntfy[count] = 0
      count += 1
    end
    return tau, duration, lmd, p, c, c_ntfy, count
  end

  def findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, f_array)
    getResultObj = Java::ReturningValues2.new
    resultObject = getResultObj.getResultFromMatlab(tau, duration, lmd, p, c, c_ntfy, f_array)
    resultObject.each do |answer|
      print answer.to_s + "\n"
    end
  end


  def myAssignment1(clause, helpers, numberOfSensors)
    @Solution = Set.new
    @Total_cost =0;
    myHash = Hash.new
    myHelperHash = Hash.new

    for j in 1..numberOfSensors
      myHash[j] =99999
    end

    for i in 0..helpers.count-1
      sensors = helpers[i].sensorList
      sensors.each do |s|
        if (myHash[s.sensorId] > s.cost)
          prev =0;
          if myHash[s.sensorId] !=99999
            prev = myHash[s.sensorId]
          end
          myHash[s.sensorId] = s.cost
          myHelperHash[s.sensorId] = i
          @Total_cost = @Total_cost + s.cost - prev;
        end
      end
    end

    temp = Hash.new

    for i in 0..helpers.count-1
      temp[i] =Set.new
    end

    myHelperHash.each do |help|
      s = temp[help[1]]
      temp[help[1]] = s<<help[0]
    end

    temp.each do |t|
      t1=t[0];
      t2 = t[1];
      if (t2.length!=0)
        @Solution.add(t2)
      end
    end

    countClause = 0
    @Solution.each do |c|
      #print c.to_a, " cost: ", @CostHash[c.to_set], "  "
      countClause +=1
      #@Total_cost += @CostHash[c.to_set]
    end

    width = 31
    f_array = Array.new(countClause) {Array.new(width)}

    for i in 0..(countClause - 1)
      for j in 0..30
        f_array[i][j] = 0
      end
    end

    index = 0

    countRow = -1
    countColumn = 0
    @Solution.each do |c|
      countRow += 1
      f_array[countRow][countColumn] = 1
      countColumn += 1
      numberOfAtoms = c.to_a.length
      f_array[countRow][countColumn] = numberOfAtoms
      countColumn += 1
      c.each do |atom|
        f_array[countRow][countColumn] = atom
        countColumn += 1
      end
      countColumn = 0
    end

    print "\n Total cost", @Total_cost
    print "\n Total number of clauses are", countClause;
    print "\n"
    f_array.each do |subarray|
      subarray.each do |element|
        print element.to_s + " "
      end
      puts " "
    end
    return f_array
  end

  #Starting point for best local search
  def makeAssignment(expression, helperList, bit)
    subClauseList = []
    possibleHelperCombinations = []

    expression = expression.flatten
    expression = expression.uniq

    #puts "EXP #{expression}"
    return findBestAssignment(expression, helperList, bit)
  end

  def makeAssignment_6(expression, helperList)
    subClauseList = []
    possibleHelperCombinations = []

    expression = expression.flatten
    expression = expression.uniq

    return findBestAssignment_6(expression, helperList)
  end

=begin
This method implements the better than greedy algorithm. Its parameters are the CNF equation, the helpers set and the k value.
It outputs the time and cost of the solution calculated by this algorithm.
=end
  def betterThan_Greedy(expression, helpers, k)

    subClauseList = []
    possibleHelperCombinations = []

    @CNFElements = expression.flatten
    @CNFElements = @CNFElements.uniq

    @K = k
    @Universe = @CNFElements.to_set
    @Subset = Set.new
    hashSet = Hash.new
    @Cost = []
    helpers.each do |c|
      @Subset << c.sensors.to_set
      hashSet[c.sensors.to_set] = c.cost
      @Cost << c.cost
    end

    @Solution = Set.new
    time = Time.new
    starttime = Time.now()
    localSearchTime = Benchmark.realtime {
      #counter = 0
=begin
    for i in @Subset do
      localcost = hashSet[i.to_set]
        set = Set.new
        for l in 1..(i.length) do
         set += i.to_a.combination(l).to_set
        end
        set.delete('[]')
        for j in set do

            if hashSet.has_key?(j.to_set)
                value = hashSet.fetch(j.to_set)
                if value > localcost
                  hashSet[j.to_set] = localcost
                end
            else
                hashSet[j.to_set] = localcost
            end
         end
        # counter = counter + 1
    end
=end
      hashSet.delete("")
      while !@Universe.empty?()

        weightJ = Hash.new
        ratioJ = Hash.new

        subsetS = hashSet.keys.to_set
        for j in subsetS
          tempSet = j.to_set & @Universe
          weightJ[j] = tempSet.length
          if tempSet.length != 0
            ratioJ[j] = (hashSet[j].to_f / tempSet.length.to_f)
          else
            ratioJ[j] = 100000
          end
        end

        keyset = hashSet.keys.to_set
        subcollectionC = Set.new
        setSj = Set.new
        setC = Set.new
        minC = 10000


        if (@Solution.size() > 0)
          for i in 2..@K do
            subcollectionC = keyset.to_a.combination(i).to_set

            for subcollection in subcollectionC

              tempSet = subcollection.to_set
              unionSet = Set.new
              costSum_subcollection = 0.0
              for j in tempSet
                #costSum_subcollection += hashSet[j]
                unionSet = unionSet | j
              end

              abc = unionSet.to_a
              for set in @Solution
                if set.to_set.subset? unionSet
                  tempSetX = unionSet & @Universe

                  if tempSetX.length > 0

                    for j in tempSet
                      costSum_subcollection += hashSet[j]
                      #unionSet = unionSet | j
                    end
                    ratio_sub = (costSum_subcollection - (hashSet[set.to_set])) / tempSetX.length
                    if ratio_sub < minC
                      minC = ratio_sub
                      setC = subcollection
                      setSj = set
                    end

                  end
                end
              end
            end
          end
        end
        #minimize between Rj and Rj,C
        tempSetY = ratioJ.keys().to_a
        setJ = tempSetY[0]
        minJ = ratioJ[setJ]

        for i in tempSetY
          if ratioJ[i] < minJ
            minJ = ratioJ[i]
            setJ = i
          elsif ratioJ[i] == minJ
            if i.length > setJ.length
              minJ = ratioJ[i]
              setJ = i
            end
          end

        end

        # greedy set
        alphaK = 1.0 - (1.0/(@K*@K*@K))

        if minJ <= minC/alphaK
          @Universe = @Universe - setJ
          @Solution << setJ
          #print setJ.to_a , ", " , minJ, "\n"
        else
          for j in setC
            @Universe -= j
          end
          #@Universe = @Universe - setC
          @Solution.delete(setSj)
          @Solution.merge(setC)
          #print "delete ", setSj.to_a , ", added " , setC , "\n"
        end


      end
    }
    endtime = Time.now()
    #print setJ.to_a , ", " , minJ
    print "\n\nBetter than greedy:  "
    total_cost = 0.0
    countClause = 0
    @Solution.each do |c|
      print c.to_a, " cost: ", hashSet[c.to_set], "  "
      total_cost += hashSet[c.to_set]
      countClause += 1
    end
    sizeOfF = countClause * 31


    width = 31
    f_array = Array.new(countClause) {Array.new(width)}

    for i in 0..(countClause - 1)
      for j in 0..30
        f_array[i][j] = 0
      end
    end

    index = 0

    countRow = -1
    countColumn = 0
    @Solution.each do |c|
      countRow += 1
      f_array[countRow][countColumn] = 1
      countColumn += 1
      numberOfAtoms = c.to_a.length
      f_array[countRow][countColumn] = numberOfAtoms
      countColumn += 1
      c.each do |atom|
        f_array[countRow][countColumn] = atom
        countColumn += 1
      end
      countColumn = 0
    end

    f_array.each do |subarray|
      subarray.each do |element|
        print element.to_s + " "
      end
      puts " "
    end

    print "\n\nTotal cost: ", total_cost.to_s
    #print "\nTotal Time: ", endtime-starttime
    print "\nTIME FOR Better than greedy Algo SEARCH #{localSearchTime}"
    return f_array;
  end

  def findBestAssignment(subClause, possibleHelperCombinations, bit)
    #find initial solution
    out = nil;
    localSearchTime = Benchmark.realtime {
      if bit == 1
        print('BLS');
        currentSolution = @@initialSolutionBLS
        initial_cost = @@initialCostBLS
      else
        print('Greedy');
        currentSolution = @@initialSolutionGreedy
        initial_cost = @@initialCostGreedy;
      end
      print "\n\nBest Local Search\n"
      puts "INITIAL COST Best Local Search #{initial_cost}"
      solution = nil
      loop do
        solution = searchForSolution(subClause, possibleHelperCombinations, currentSolution)
        if solution == nil
          print("No neighbor found!!\n");
          solution = currentSolution
          break;
        end
        final_arr = []
        final_cost = 0
        solution.each do |current_helper|
          final_cost = final_cost + current_helper.cost
          final_arr << current_helper.sensors
        end
        break if (initial_cost <= final_cost)
        currentSolution = solution[0..solution.size-1]
        initial_cost = final_cost
      end
      final_arr = []
      final_cost = 0
      countClause = 0
      countAtoms = 0
      solution.each do |current_helper|
        final_cost = final_cost + current_helper.cost
        final_arr << current_helper.sensors
        countClause += 1
      end

      puts "number of clauses = " + countClause.to_s

      sizeOfF = countClause * 31

      width = 31
      f_array = Array.new(countClause) {Array.new(width)}

      for i in 0..(countClause - 1)
        for j in 0..30
          f_array[i][j] = 0
        end
      end

      index = 0

      countRow = -1
      countColumn = 0
      solution.each do |c|
        countRow += 1
        f_array[countRow][countColumn] = 1
        countColumn += 1
        numberOfAtoms = c.sensors.to_a.length
        f_array[countRow][countColumn] = numberOfAtoms
        countColumn += 1
        c.sensors.each do |atom|
          f_array[countRow][countColumn] = atom
          countColumn += 1
        end
        countColumn = 0
      end

      f_array.each do |subarray|
        subarray.each do |element|
          print element.to_s + " "
        end
        puts " "
      end
      out = f_array;
      puts "FINAL COST for Best Local Search #{final_cost}"
    }
    puts "TIME FOR Best Local Search SEARCH #{localSearchTime}"
    return out;
  end

  def searchForSolution(expression, possibleHelperCombinations, currentSolution)
    helperListCopy = Array.new
    currentSolutionCopy = Array.new
    currentSolutionCopy = currentSolution[0..currentSolution.size-1]
    helperListCopy = possibleHelperCombinations[0..possibleHelperCombinations.size-1]
    neighbors = findNeighbors(expression, helperListCopy, currentSolutionCopy)
    bestNeighbor = findBestNeighbor(neighbors)
    return bestNeighbor
  end

  def findBestNeighbor(neighbors)
    lowestCost = 999999999999
    bestNeighbor = nil
    neighbors.each do |neighbor|
      neighborCost = 0
      neighbor.each do |helper|
        neighborCost+= helper.cost
      end
      if (neighborCost < lowestCost)
        bestNeighbor = neighbor
        lowestCost = neighborCost
      end
    end
    return bestNeighbor
  end

  def findNeighbors(expression, helperListCopy, currentSolutionCopy)
    finalNeighbors = Array.new
    tempCurrentSolution = Array.new
    i = 0
    currentSolutionCopy.each do |solHelper|
      helperListCopy.each do |helper|
        tempCurrentSolution = currentSolutionCopy[0..currentSolutionCopy.size-1]
        tempCurrentSolution.delete(solHelper)
        tempCurrentSolution << helper
        sensors = []
        tempCurrentSolution.each do |h|
          sensors << h.sensors
        end
        sensors = sensors.flatten
        if (expression-sensors).empty?
          finalNeighbors << tempCurrentSolution
        end
      end
    end
    return finalNeighbors
  end


  def findInitialSolution(subClause, possibleHelperCombinations)
    subClause = subClause.flatten
    subClause = subClause.uniq
    currentSolution = []
    loop do
      #GIVE NUMBER OF HELPERS TO CONSIDER FOR FIRST RANDOM ASSIGNMENT
      randomStart = Random.rand(0..possibleHelperCombinations.size-1); #need to change this TODO
      lowerBound = possibleHelperCombinations.size/4
      upperBound = possibleHelperCombinations.size/3
      randomStart = Random.rand(lowerBound..upperBound)
      #randomStart = Random.rand(10..15)
      sensors = []
      selectedHelperCombination = possibleHelperCombinations.sample(randomStart)
      selectedHelperCombination.each do |helper|
        sensors << helper.sensors
      end
      sensors = sensors.flatten
      #check if all the required sensors are covered by the selected helpers
      if (subClause-sensors).empty?
        currentSolution = selectedHelperCombination;
        #puts "RandomStartSeed = #{randomStart}"
        return currentSolution

      end
    end
  end

=begin
This method implements the Tabu search heuristic.
=end
  #expression = @@clause, helperList = @@helpers
  def makeAssignment_5(expression, helperList, bit)
    subClauseList = []
    possibleHelperCombinations = []

    expression = expression.flatten
    expression = expression.uniq

    #puts "EXP #{expression}"
    return findBestAssignment_5(expression, helperList, bit)
  end

  def findBestAssignment_6(subClause, possibleHelperCombinations)
    #find initial solution
    out = nil;
    localSearchTime = Benchmark.realtime {
      count = 0
      fInitialCost = 0
      fCost = 999999999999999999
      fSolution = []
      countClause = 0
      countAtoms = 0
      finalsolution = nil
      loop do
        currentSolution = findInitialSolution_ls(subClause, possibleHelperCombinations)
        initial_cost = 0
        current_arr = []
        currentSolution.each do |current_helper|
          initial_cost = initial_cost + current_helper.cost
          current_arr << current_helper.sensors
        end
        if fInitialCost == 0
          fInitialCost = initial_cost
        end
        solution = nil
        loop do
          solution = searchForSolution(subClause, possibleHelperCombinations, currentSolution)
          final_arr = []
          final_cost = 0
          solution.each do |current_helper|
            final_cost = final_cost + current_helper.cost
            final_arr << current_helper.sensors
          end
          break if (initial_cost <= final_cost)
          currentSolution = solution[0..solution.size-1]
          initial_cost = final_cost
        end
        final_arr = []
        final_cost = 0
        countClause = 0
        solution.each do |current_helper|
          final_cost = final_cost + current_helper.cost
          final_arr << current_helper.sensors

        end
        if final_cost >= fCost
          count +=1
        end
        if final_cost < fCost
          fCost = final_cost
          fsolution = final_arr
        end

        if count == 20
          finalsolution = solution
          break
        end
      end

      finalsolution.each do |c|
        countClause += 1
      end

      puts "number of clause = " + countClause.to_s

      sizeOfF = countClause * 31


      width = 31
      f_array = Array.new(countClause) {Array.new(width)}

      for i in 0..(countClause - 1)
        for j in 0..30
          f_array[i][j] = 0
        end
      end

      index = 0

      countRow = -1
      countColumn = 0
      finalsolution.each do |c|
        countRow += 1
        f_array[countRow][countColumn] = 1
        countColumn += 1
        numberOfAtoms = c.sensors.to_a.length
        f_array[countRow][countColumn] = numberOfAtoms
        countColumn += 1
        c.sensors.each do |atom|
          f_array[countRow][countColumn] = atom
          countColumn += 1
        end
        countColumn = 0
      end

      f_array.each do |subarray|
        subarray.each do |element|
          print element.to_s + " "
        end
        puts " "
      end
      out = f_array;

      print "\n\nMultiple Start Best Local Search"
      puts "\nINITIAL COST for Multiple start BLS #{fInitialCost}"
      print "\nFINAL COST for Multiple start BLS  #{fCost}"
    }
    puts "\nTIME FOR Multiple start best local SEARCH #{localSearchTime}"
    return out;
  end

=begin
This method implements the Tabu search heuristic.
=end
  def findBestAssignment_5(subClause, possibleHelperCombinations, bit)
    #find initial solution
    out = nil
    fifthSearchTime = Benchmark.realtime {
      if bit == 1
        print('NORMAL');
        currentSolution = findInitialSolution_ls(subClause, possibleHelperCombinations)
      else
        print('Greedy');
        currentSolution = @@initialSolutionGreedy
        initial_cost = @@initialCostGreedy;
      end

      #currentSolution = initialSolution
      #TODO
      currentSolution.sort! {|a, b| a.id <=> b.id}
      #puts "INITIAL SOLUTION #{currentSolution}"

      @previousSolutionHash = Hash.new
      currentSolString = String.new
      #currentSolution.each do |a|
      #  currentSolString = currentSolString + "-" + a.id.to_s
      #end
      @previousSolutionHash[currentSolString] = true
      # @previousSolutions = Array.new
      # @previousSolutions << currentSolution
      #puts "INITIAL SOLUTION #{currentSolution}"
      initial_cost = 0
      current_arr = []
      currentSolution.each do |current_helper|
        initial_cost = initial_cost + current_helper.cost
        current_arr << current_helper.sensors
      end
      finalBestSolution = currentSolution
      finalBestCost = initial_cost
      tempInitialCost = initial_cost
      #puts "INITIAL COST #{initial_cost}"
      if (subClause - current_arr.flatten).empty?
        #puts "TRUE INITIAL"
      else
        #puts "FALSE INITIAL"
      end
      solution = nil
      loopCount = 0
      finalsolution = nil
      loop do
        loopCount += 1
        solution = searchForSolution_5(subClause, possibleHelperCombinations, currentSolution)
        final_arr = []
        final_cost = 0

        tempCost = 99999999

        solution.each do |current_helper, cost|
          if (cost < tempCost)
            tempCost = cost
            currentSolution = current_helper
          end
        end

        if (tempCost < finalBestCost)

          finalBestSolution = currentSolution
          finalBestCost = tempCost
        end


        currentSolution.sort! {|a, b| a.id <=> b.id}
        currentSolString = String.new
        currentSolution.each do |a|
          currentSolString = currentSolString + "-" + a.id.to_s
        end
        @previousSolutionHash[currentSolString] = true
        #puts("COSTS", initial_cost, final_cost)
        if (loopCount == (@@helpers.size))
          puts finalBestSolution
          finalsolution = finalBestSolution
          break
        end
        #currentSolution = solution[0..solution.size-1]
        #initial_cost = final_cost
      end

      countClause = 0

      finalsolution.each do |helper|
        countClause += 1
      end

      puts "number of clasues = " + countClause.to_s

      width = 31
      f_array = Array.new(countClause) {Array.new(width)}

      for i in 0..(countClause - 1)
        for j in 0..30
          f_array[i][j] = 0
        end
      end

      index = 0

      countRow = -1
      countColumn = 0
      finalsolution.each do |c|
        countRow += 1
        f_array[countRow][countColumn] = 1
        countColumn += 1
        numberOfAtoms = c.sensors.to_a.length
        f_array[countRow][countColumn] = numberOfAtoms
        countColumn += 1
        c.sensors.each do |atom|
          f_array[countRow][countColumn] = atom
          countColumn += 1
        end
        countColumn = 0
      end

      f_array.each do |subarray|
        subarray.each do |element|
          print element.to_s + " "
        end
        puts " "
      end
      out = f_array;
      print "\n5th algoritm Search results \n"
      puts "Initial COST #{tempInitialCost}"
      puts "Final COST #{finalBestCost}"
      #puts("PREVIOUS SIZE", @previousSolutionHash.length)
    }
    puts "TIME FOR 5th SEARCH #{fifthSearchTime}\n"
    return out;
  end


=begin
This method implements the Tabu search heuristic  in finding its neighbours and the best amoung them.
=end

  def searchForSolution_5(expression, possibleHelperCombinations, currentSolution)
    helperListCopy = Array.new
    currentSolutionCopy = Array.new
    #puts("SDFGSDFGSDFGSDFGDSFG", currentSolution.inspect)
    currentSolutionCopy = currentSolution[0..currentSolution.size-1]
    helperListCopy = possibleHelperCombinations[0..possibleHelperCombinations.size-1]
    #puts("HELPER LIST", helperListCopy.inspect)
    neighbors = findNeighbors_5(expression, helperListCopy, currentSolutionCopy)
    bestNeighbor = findBestNeighbor_5(neighbors)
    #puts("BEST NEIGHBOR", bestNeighbor)
    return bestNeighbor
  end

  def findBestNeighbor_5(neighbors)
    #puts("NEIGFHBOT SIZER", neighbors.size)
    lowestCost = 999999999999
    bestNeighbors = Hash.new
    neighbors.each do |neighbor|
      neighborCost = 0
      neighbor.each do |helper|
        neighborCost+= helper.cost
      end
      bestNeighbors[neighbor] = neighborCost
    end
    bestNeighbors = bestNeighbors.sort_by {|x, y| y}
    #puts("BEST NEIGHBORS", bestNeighbors.take(2))
    return bestNeighbors
  end

  #TODO
  def isNewSolution(currentSolutionTemp)
    currentSolutionTemp.sort! {|a, b| a.id <=> b.id}
    currentSolString = String.new
    currentSolution.each do |a|
      currentSolString = currentSolString + "-" + a.id.to_s
    end
    return @previousSolutionHash.has_key?(currentSolString)

  end

  def findNeighbors_5(expression, helperListCopy, currentSolutionCopy)
    finalNeighbors = Array.new
    tempCurrentSolution = Array.new
    i = 0
    currentSolutionCopy.each do |solHelper|
      helperListCopy.each do |helper|
        tempCurrentSolution = currentSolutionCopy[0..currentSolutionCopy.size-1]
        tempCurrentSolution.delete(solHelper)
        tempCurrentSolution << helper
        sensors = []
        tempCurrentSolution.each do |h|
          sensors << h.sensors
        end
        sensors = sensors.flatten

        tempCurrentSolution.sort! {|a, b| a.id <=> b.id}
        currentSolString = String.new
        tempCurrentSolution.each do |a|
          currentSolString = currentSolString + "-" + a.id.to_s
        end

        #TODO check both if conditions to compare what works
        #if((expression-sensors).empty? && !(@previousSolutions.include?(tempCurrentSolution)))
        #if((expression-sensors).empty? && isNewSolution(tempCurrentSolution))
        if ((expression-sensors).empty? && !(@previousSolutionHash.has_key?(currentSolString)))
          finalNeighbors << tempCurrentSolution
        end
      end
    end
    #puts("FINAL NEIGHBORS", finalNeighbors.inspect)
    #puts("SIZE", finalNeighbors.size)
    return finalNeighbors
  end

  def findInitialSolution_5(subClause, possibleHelperCombinations)
    currentSolution = []
    loop do
      #GIVE NUMBER OF HELPERS TO CONSIDER FOR FIRST RANDOM ASSIGNMENT
      randomStart = Random.rand(0..possibleHelperCombinations.size-1); #need to change this TODO
      lowerBound = possibleHelperCombinations.size/4
      upperBound = possibleHelperCombinations.size/3
      randomStart = Random.rand(lowerBound..upperBound)
      #randomStart = Random.rand(10..15)
      sensors = []
      selectedHelperCombination = possibleHelperCombinations.sample(randomStart)
      selectedHelperCombination.each do |helper|
        sensors << helper.sensors
      end
      sensors = sensors.flatten
      #check if all the required sensors are covered by the selected helpers
      if (subClause-sensors).empty?
        currentSolution = selectedHelperCombination;
        #puts "RandomStartSeed = #{randomStart}"
        return currentSolution
      end
    end
  end

=begin
  This method implements the greedy set cover algorithm and its
  parameters are the CNF clause and the helpers set.
  o/p: Time and cost of the solution computed by greedy algorithm.
=end

  def greedySetCover(clause, helpers)

    subClauseList = []
    possibleHelperCombinations = []
    @CNFClause = clause.to_a
    @CNFElements = []
    @CNFElements = clause.flatten
    @CNFElements = @CNFElements.uniq.to_set

    @Universe = @CNFElements.to_set

    @Total_cost = 0
    @Collapse = Set.new
    @Solution = Set.new
    @CostHash = Hash.new
    @Subset = Set.new
    @Cost = []
    @tempArray = Array.new
    @HelperMap = Hash.new

    helpers.each do |c|
      @Subset << c.sensors.to_set
      @CostHash[c.sensors.to_set] = c.cost
      @Cost << c.cost
      @HelperMap[c.sensors.to_set] = c;
    end


    starttime = Time.now()
    counter = 0


    temp_included = Hash.new
    while !solution_size

      len, j = 100000, 10000
      for subset in @Subset do

        j = check_elements(subset.to_a)
        if !temp_included.has_key?(subset) && j < len # && set_len < temp.length
          index = subset #set_len = i, temp.length
          len = j
        end

      end
      temp_included[index] = true
      tempCollapse = Array.new
      t= Set.new;
      index.each do |a|
        t.add(a)
      end

      @Sol = Set.new
      @Solution.each do |b|
        @Sol += b
      end
      @Total_cost = @Total_cost + @CostHash[index.to_set]
      index = index.to_set

      @Sol.each do |b|
        index.delete(b)
      end
      @CNFClause.each do |set|
        if (index.to_set.intersection(set.to_set))
          tempCollapse << (index & set).to_a
        end
      end
      #tempCollapse << (index - tempCollapse).to_a
      @Solution.add(index)
      if (@HelperMap[index]!=nil)
        @tempArray << @HelperMap[index];
      end
      @Collapse << tempCollapse.to_set
    end
    @@initialSolutionGreedy = @tempArray[0..@tempArray.size-1]
    @@initialCostGreedy = @Total_cost

    print "\n\nJust greedy:  "
    endtime = Time.now()
    countClause = 0
    @Solution.each do |c|
      print c.to_a, " cost: ", @CostHash[c.to_set], "  "
      countClause +=1
      #@Total_cost += @CostHash[c.to_set]
    end
    print "\n collapsed output: ", @Collapse.to_a
    print "\n CNF:  ", clause.to_a
    print "\nTotal cost: ", @Total_cost.to_s
    print "\nTotal Time: ", endtime-starttime
    puts "number of clause " + countClause.to_s

    sizeOfF = countClause * 31

    width = 31
    f_array = Array.new(countClause) {Array.new(width)}

    for i in 0..(countClause - 1)
      for j in 0..30
        f_array[i][j] = 0
      end
    end

    index = 0

    countRow = -1
    countColumn = 0
    @Solution.each do |c|
      countRow += 1
      f_array[countRow][countColumn] = 1
      countColumn += 1
      numberOfAtoms = c.to_a.length
      f_array[countRow][countColumn] = numberOfAtoms
      countColumn += 1
      c.each do |atom|
        f_array[countRow][countColumn] = atom
        countColumn += 1
      end
      countColumn = 0
    end

    f_array.each do |subarray|
      subarray.each do |element|
        print element.to_s + " "
      end
      puts " "
    end
    return f_array;
  end

=begin
This method calculates the total number of newly added elements in every step of greedy algorithm.
=end
  def check_elements(subset)
    temp = Set.new
    @Sol = Set.new
    @Solution.each do |b|
      @Sol += b
    end

    temp = subset.to_set

    @Sol.each do |b|
      temp.delete(b)
    end
    count = 0.0
    count = @CostHash[subset.to_set]

    if temp.length == 0
      return 100000
    else
      count = count.to_f / temp.size()
      return count
    end
  end

=begin
This method checks if a solution if found.
=end
  def solution_size
    temp = Set.new
    @Solution.each do |c|
      temp = temp | c
    end
    return (@CNFElements - temp).empty?()

  end

=begin
 Local Search
=end

  def makeAssignment_ls(expression, helperList)
    subClauseList = []
    possibleHelperCombinations = []

    #Determine unique sensors in expression
    expression = expression.flatten
    expression = expression.uniq

    #puts "EXP #{expression}"
    return findBestAssignment_ls(expression, helperList)
  end

  def findBestAssignment_ls(subClause, possibleHelperCombinations)
    #find initial solution
    out = nil;
    localSearchTime = Benchmark.realtime {
      currentSolution = findInitialSolution_ls(subClause, possibleHelperCombinations)
      initial_cost = 0
      current_arr = []
      currentSolution.each do |current_helper|
        initial_cost = initial_cost + current_helper.cost
        current_arr << current_helper.sensors
      end
      @@initialSolutionBLS = currentSolution[0..currentSolution.size-1]
      @@initialCostBLS = initial_cost
      solution = searchForSolutionUsingSA(subClause, currentSolution)
      solution = solution.uniq
      #puts "Best helper assignment for SUBCLAUSE using local search #{subClause} = #{solution}"
      final_arr = []
      final_cost =
          countClause = 0
      countAtoms = 0
      solution.each do |current_helper|
        final_cost = final_cost + current_helper.cost
        final_arr << current_helper.sensors
        countClause += 1
      end
      puts "number of clauses = " + countClause.to_s

      sizeOfF = countClause * 31


      width = 31
      f_array = Array.new(countClause) {Array.new(width)}

      for i in 0..(countClause - 1)
        for j in 0..30
          f_array[i][j] = 0
        end
      end

      index = 0

      countRow = -1
      countColumn = 0
      solution.each do |c|
        countRow += 1
        f_array[countRow][countColumn] = 1
        countColumn += 1
        numberOfAtoms = c.sensors.to_a.length
        f_array[countRow][countColumn] = numberOfAtoms
        countColumn += 1
        c.sensors.each do |atom|
          f_array[countRow][countColumn] = atom
          countColumn += 1
        end
        countColumn = 0
      end

      f_array.each do |subarray|
        subarray.each do |element|
          print element.to_s + " "
        end
        puts " "
      end
      out = f_array;
      print "\nSimulated Annealing (Local Search)\n"
      print "\nINITIAL COST of Simulated Annealing (Local Search) #{initial_cost}"
      print "\nFINAL COST #{final_cost}"
    }
    print "\nTIME FOR Simulated Annealing (LOCAL SEARCH) #{localSearchTime}"
    return out;

  end

  def searchForSolutionUsingSA(expression, currentSolution)
    tempSol = Array.new
    localSearchTime = Benchmark.realtime {
      expression = expression.flatten
      expression = expression.uniq
      tempSol = currentSolution
      helperSize = @@helpers.size
      kMax = (helperSize*(helperSize-1))/2
      (kMax-1).downto(1) do |k|
        temperature = k.to_f/kMax
        #puts("TEMPERATURE", temperature)
        randomNeighbor = findRandomNeighborSA(expression, tempSol[0..tempSol.size-1])
        solutionCost = 0
        neighborCost = 0
        randomNeighbor.each do |neighbor|
          neighborCost += neighbor.cost
        end
        tempSol.each do |fSolution|
          solutionCost += fSolution.cost
        end
        if (neighborCost < solutionCost)
          tempSol = randomNeighbor
        else
          probability = Math.exp(-((neighborCost-solutionCost)/temperature))
          if probability > Random.rand() #(Random.rand(0..100) < 10)#
            tempSol = randomNeighbor
          end
        end
      end
    }
    #puts "Time taken for local Search = #{localSearchTime}"
    return tempSol
  end

  def findRandomNeighborSA(expression, runninSol)
    tempSolution = Array.new
    runninSolSize = runninSol.size
    tempSolution[0..runninSolSize] = runninSol[0..runninSolSize]
    loopCount = tempSolution.size*@@helpers.size
    for i in 0..loopCount
      randSolution = Random.rand(0..tempSolution.size-1)
      randHelper = Random.rand(0..@@helpers.size-1)
      tempSolution[randSolution] = @@helpers[randHelper]
      temp = []
      tempSolution.each do |sol|
        temp << sol.sensors
      end
      if (expression-temp.flatten).empty?
        return tempSolution
      else
        tempSolution[0..runninSolSize] = runninSol[0..runninSolSize]
      end
    end
    return tempSolution
  end

  def findInitialSolution_ls(subClause, possibleHelperCombinations)
    currentSolution = []
    lowerBound = 1
    # start from 1 helper and check if all sensors are covered, if not, increase helper size and keep checking till
    #solution is found
    loop do
      #GIVE NUMBER OF HELPERS TO CONSIDER FOR FIRST RANDOM ASSIGNMENT
      upperBound = possibleHelperCombinations.size
      randomStart = lowerBound
      sensors = []
      selectedHelperCombination = possibleHelperCombinations.sample(randomStart)
      selectedHelperCombination.each do |helper|
        sensors << helper.sensors
      end
      sensors = sensors.flatten
      #check if all the required sensors are covered by the selected helpers
      if (subClause-sensors).empty?
        currentSolution = selectedHelperCombination;
        return currentSolution
      else
        lowerBound = (lowerBound+1)%upperBound
      end
    end
  end


  def makeAssignment_blsVariation(subClause, helperList)
    blsVariationTime = Benchmark.realtime {
      currentSolution = @@initialSolutionBLS[0..@@initialSolutionBLS.size-1]
      initialCost = @@initialCostBLS

      #FIND BEST POSSIBLE SOLUTION
      solution = searchForSolutionUsingBLSVariation(subClause, currentSolution)


      solution = solution.uniq
      #puts "Best helper assignment for SUBCLAUSE using local search #{subClause} = #{solution}"
      final_arr = []
      final_cost = 0
      solution.each do |current_helper|
        final_cost = final_cost + current_helper.cost
        final_arr << current_helper.sensors
      end

      if (subClause - final_arr.flatten).empty?
        #puts "TRUE"
      else
        #puts "FALSE"
      end

      print "\nBEST LOCAL SEARCH VARIATION\n"
      print "\nINITIAL COST of BEST Local Search VARIATION #{initialCost}"
      print "\nFINAL COST #{final_cost}"

    }
    puts("\nTime for variation of best local search", blsVariationTime)
  end

  def searchForSolutionUsingBLSVariation(expression, currentSolution)
    tempSol = Array.new
    localSearchTime = Benchmark.realtime {
      expression = expression.flatten
      expression = expression.uniq
      tempSol = currentSolution[0..currentSolution.size-1]
      helperSize = @@helpers.size
      kMax = (helperSize*(helperSize-1))/2
      for k in 1..kMax-1
        temperature = k.to_f/kMax

        randomNeighbor = findRandomNeighborBLSVariation(expression, tempSol[0..tempSol.size-1])
        solutionCost = 0
        neighborCost = 0
        randomNeighbor.each do |neighbor|
          neighborCost += neighbor.cost
        end
        tempSol.each do |fSolution|
          solutionCost += fSolution.cost
        end
        if (neighborCost < solutionCost)
          tempSol = randomNeighbor
        else
          probability = Math.exp(-((neighborCost-solutionCost)/temperature))
          if probability > Random.rand() #(Random.rand(0..100) < 10)#
            tempSol = randomNeighbor
          end
        end
      end
    }
    #puts "Time taken for local Search = #{localSearchTime}"
    return tempSol
  end

  def findRandomNeighborBLSVariation(expression, currentSolution)
    tempSolution = Array.new
    currentSolutionSize = currentSolution.size
    tempSolution[0..currentSolutionSize] = currentSolution[0..currentSolutionSize]
    worstSensor = Sensor.new(-1, -1) #keep track of worst sensor
    worstHelperID = -1 #keep track of helper with worst sensor

    #FIND WORST SENSOR IN THE CURRENT SOLUTION  (First, select a sensor with worst cost)
    currentSolution.each do |curSolHelper|
      sensorList = curSolHelper.sensorList
      #puts("SENSOR LIST", sensorList.inspect)
      sensorList.each do |sensor|
        if (sensor.cost > worstSensor.cost)
          worstSensor = sensor
          worstHelperID = curSolHelper.id
        end
      end
    end

    #puts("WORST SENSOR", worstSensor.inspect)
    #puts("WORSR HELPER ID", worstHelperID)

    #Find the best helper with smallest cost for the sensor
    bestSensor = Sensor.new(-1, 99999)
    bestHelper = nil
    tempHelperList = @@helpers[0..@@helpers.size-1]
    tempHelperList.each do |indHelper|
      #puts("IND HELPER", indHelper.inspect)
      if (indHelper.sensors.include?(worstSensor.sensorId))
        #puts("IN TRUE 1")
        tempSensor = indHelper.sensorList.find {|s| s.sensorId == worstSensor.sensorId}
        #puts("TEMP SENSOR ", tempSensor.inspect)
        if ((tempSensor.cost < worstSensor.cost) && (tempSensor.cost < bestSensor.cost) && areAllSensorsCovered(expression, indHelper, currentSolution, worstHelperID))
          #puts("IN TRUE 2")
          bestSensor = tempSensor
          bestHelper = indHelper
        end
      end
    end

    #REPLACE WORST SENSOR-HELPER WITH BEST SENSOR-HELPER
    if (bestHelper != nil)
      tempSolution.delete_if {|s| s.id == worstHelperID}
      tempSolution << bestHelper
    end
    return tempSolution
=begin
    loopCount = tempSolution.size*@@helpers.size
    for i in 0..loopCount
      randSolution = Random.rand(0..tempSolution.size-1)
      randHelper = Random.rand(0..@@helpers.size-1)
      tempSolution[randSolution] = @@helpers[randHelper]
      temp = []
      tempSolution.each do |sol|
        temp << sol.sensors
      end
      if (expression-temp.flatten).empty?
        return tempSolution
      else
        tempSolution[0..runninSolSize] = runninSol[0..runninSolSize]
      end
    end
    return tempSolution
=end
  end

  def areAllSensorsCovered(expression, indHelper, currentSolution, worstHelperID)
    tempSensorList = Array.new
    tempSensorList << indHelper.sensorList
    tempCurrentSolution = currentSolution[0..currentSolution.size-1]
    tempCurrentSolution.delete_if {|s| s.id == worstHelperID}
    tempCurrentSolution.each do |t|
      #puts("TTTTT", t.inspect)
      tempSensorList << t.sensors
    end
    #puts("EXPRESSION", expression.inspect)
    #puts("TEMP SENSOR LIST", tempSensorList.inspect)
    if (expression - tempSensorList.flatten).empty?
      return true
    else
      return false
    end

  end


end

=begin
 Each sensor object has a sensorId and a cost
=end
class Sensor
  attr_accessor :cost, :sensorId

  def initialize(sensorId, cost)
    @sensorId = sensorId
    @cost = cost
  end
end

=begin
 Each Helper has a list of sensors, total cost and ID
=end
class Helper
  attr_accessor :sensorList, :cost, :sensors, :id

  def initialize(sensorList, sensors, cost, id)
    @sensors = sensors
    @cost = cost
    @id = id
    @sensorList = sensorList
  end
end

helpers = Array.new;
#set the number of helpers and sensors for tests here and the range of costs for sensors
numberOfHelpers = 90
numberOfSensors = 185
sensorCostMin = 50
sensorCostMax = 120
randomSensors = Random.new
listOfSensors = *(1..numberOfSensors)
listOfSensors = listOfSensors.shuffle
tempHelper = Array.new
individualSensorCost = Array.new # Array to store individual cost of sensors
sensorList = Array.new
helperCost = 0
for i in 0..numberOfHelpers # NO OF HELPERS
  numberOfSensorsInHelper = Random.rand(3..8) # NO OF SENSORS IN HELPER
  tempHelper.push(listOfSensors.pop(numberOfSensorsInHelper))
  tempHelper = tempHelper.flatten
  for j in 0..numberOfSensorsInHelper-1
    temp = Random.rand(sensorCostMin..sensorCostMax)
    sensorList << Sensor.new(tempHelper[j], temp)
    individualSensorCost << temp
  end
  sensorList.delete_if {|s| s.sensorId == nil || s.cost == nil}
  helpers << Helper.new(sensorList, tempHelper, individualSensorCost.inject(:+), i) #CREATE HELPER WITH COST
  sensorList = []
  tempHelper = []
  if (listOfSensors.size < 9)
    listOfSensors = *(1..numberOfSensors)
    listOfSensors = listOfSensors.shuffle
  end
end

assignHelpers = Assign.create
clause = []
tempSubClause = []
randomSensors = Random.new
for i in 0..numberOfHelpers
  temp = Random.rand(3..8)
  for j in 0..temp
    tempSubClause << randomSensors = Random.new.rand(1..numberOfSensors) #NO OF SENSORS  #GENERATE SENSOR THAT GOES IN A HELPER
  end
  clause << tempSubClause
  tempSubClause = []
end
@@helpers = Array.new;
@@helpers = helpers

#puts "The clause is as follows:"
#puts clause

tau, duration, lmd, p, c, c_ntfy, count = assignHelpers.init(numberOfSensors);


print("----------FINE GREEDY------------\n")
out4 = assignHelpers.myAssignment1(clause, helpers, numberOfSensors)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out4);

#Call greedy
print("----------GREEDY ALGORITHM------------\n")
out1 = assignHelpers.greedySetCover(clause, helpers)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out1);
#assignHelpers.betterThan_Greedy(clause, helpers, 3)

#call better than greedy
#print("----------BETTER THAN GREEDY------------\n")
#out3 = assignHelpers.betterThan_Greedy(clause, helpers, 4)
#assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out3);

#Call simulated annealing (local search)
print("----------SIMULATED ANNEALING ALGORITHM------------\n")
out2= assignHelpers.makeAssignment_ls(clause, @@helpers)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out2);

#Call best local search
print("----------BEST LOCAL SEARCH - BLS ------------\n")
out5 = assignHelpers.makeAssignment(clause, helpers, 1)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out5);

#Call best local search
print("----------BEST LOCAL SEARCH - GREEDY ------------\n")
out6 = assignHelpers.makeAssignment(clause, helpers, 2)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out6);

#call tabu search
print("----------TABU ------------\n")
out7 = assignHelpers.makeAssignment_5(clause, helpers, 1)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out7);

#call tabu search
print("----------TABU GREEDY ------------\n")
out8 = assignHelpers.makeAssignment_5(clause, helpers, 2)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out8);

#Call multiple start best local search
print("----------MULTIPLE START BLS ------------\n")
out9 = assignHelpers.makeAssignment_6(clause, helpers)
assignHelpers.findMonitoringCost(tau, duration, lmd, p, c, c_ntfy, count, out9);
