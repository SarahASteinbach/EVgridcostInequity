%%Load/determine all power flow input data as saved in the workspace files,
%%no need to run to perform calculations, only included to improve
%%understanding of code

%% Select PATHs 
    pathW19 = "Y:\TUM-PC\Modelling Sarah\Daten\Mobilitätsdaten\W19.csv";
    pathH19 = "Y:\TUM-PC\Modelling Sarah\Daten\Mobilitätsdaten\HH19.csv";
    pathT19 = "Y:\TUM-PC\Modelling Sarah\Daten\Mobilitätsdaten\Tank20.csv";
    pathPV = "Y:\TUM-PC\Modelling Sarah\Daten\PV Data\";
    pathsHHLP = ["Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_December" 
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_Mar"
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_Jun"
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_December2"
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_Mar2"
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_Jun2"
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_September"
    "Y:\TUM-PC\Modelling Sarah\Daten\Household Load Profiles\Households_September2"];

%% Load mobility and PV datasets 
% Load KIT MOP Trip and household datasets
W19 = readtable(pathW19);
H19 = readtable(pathH19);
T19 = readtable(pathT19);

%% Load household load profiles
% Paths specifies the folders where all subfolders of each household are located
% order: December, March, June
paths = pathsHHLP;
for j = 1:8
    list = dir(paths(j));
    list = list([list.isdir]); 
    list = list(~ismember({list.name},{'.' '..'}));
    householdprofiles = zeros(length(list)-3,10081);
    for i=1:length(list)-3
        filepath = paths(j) + "\"+list(i).name+ "\Results\" + "Overall.SumProfiles.Electricity.csv";
        householdprofiles(i,1) = identifyHousehold(list(i).name);
        profiles = readmatrix(filepath);
        householdprofiles(i,2:end) = profiles(:,3).';
    end
    if j==1
        householdprofiles_december = householdprofiles;
    elseif j==2
        householdprofiles_march = householdprofiles;
    elseif j==3
        householdprofiles_june = householdprofiles;   
    elseif j==4
        householdprofiles_december = [householdprofiles_december; householdprofiles];
    elseif j==5
        householdprofiles_march = [householdprofiles_march; householdprofiles];
    elseif j==6
        householdprofiles_june = [householdprofiles_june; householdprofiles];        
    elseif j==7
        householdprofiles_september = householdprofiles;
    else
        householdprofiles_september = [householdprofiles_september; householdprofiles];
        clear householdprofiles   
    end
end
%% Data preprocessing of the trip dataset
% concat ID + PERSNR to create a unique identifier for each person
W19.UUID = str2double(strcat(string(int64(W19.ID)),string(int64(W19.PERSNR))));

%remove all Households where a member conducted at least one trip to hotel/vacation
W = W19;
householdHotel = W(W.ZWECK == 9,:).ID;
W = W(~ismember(W.ID,householdHotel),:);

%remove all Persons with distances > 200km and durations > 132 min
%which is 99% case
quantile(W.DAUER,0.99);
persLongDuration = W(W.DAUER > 132,:).UUID;
W = W(~ismember(W.UUID,persLongDuration),:);

persLongDistance = W(W.KM > 200,:).UUID;
W = W(~ismember(W.UUID,persLongDistance),:);

%% Split into above and below average income
H_rich=H19.ID(H19.EINKO>=8);
H_poor=H19.ID(H19.EINKO<8);
HH_rich=H19(ismember(H19.ID,H_rich),:);
HH_poor=H19(ismember(H19.ID,H_poor),:);

W_rich=W(ismember(W.ID,H_rich),:);
W_poor=W(ismember(W.ID,H_poor),:);




%% Change depending on wanting to greate the dataset for above or below average income households


%HH_group='poor';
HH_group='rich';

if HH_group=='poor'
    W=W_poor;
    H=HH_poor;
    car_segment_distribution = [0.0815 0.1917 0.2096 0.1183 0.3147 0.0842];
else
    W=W_rich;
    H=HH_rich;
    %test effect usage
    %W=W_poor;
    %H=HH_poor;

    car_segment_distribution = [0.0610 0.1340 0.1874 0.1697 0.2546 0.1933];
    
    %test effect model choice
    %car_segment_distribution = [0.0815 0.1917 0.2096 0.1183 0.3147 0.0842];
end

 

%% create a trip dataset with trips that only use car as main mode of transport
car_trips = W;
%drop all trips without car
car_trips = car_trips(car_trips.VMDIW == 41,:);
%convert a ab- and anzeit in minute indices
car_trips.abzeit_i = floor(car_trips.ABZEIT/100)*60 + mod(car_trips.ABZEIT,100);
car_trips.anzeit_i = floor(car_trips.ANZEIT/100)*60 + mod(car_trips.ANZEIT,100);
%% Create 5 minute timesteps
% round ab- and anzeit to nearest 5 min index and make sure that no trips
% are 0 minutes long
timestep_size = 5;
day_length=1440/timestep_size;
car_trips.abzeit_i_new = round(car_trips.abzeit_i/timestep_size)*timestep_size;
car_trips.anzeit_i_new = round(car_trips.anzeit_i/timestep_size)*timestep_size;
% find durations of 0
idx = find(car_trips.anzeit_i_new-car_trips.abzeit_i_new==0);
% shift half the 0 durations: start time 5 min earlier
rnd_idx = idx(randperm(numel(idx),ceil(numel(idx)/2)));
car_trips.abzeit_i_new(rnd_idx) = car_trips.abzeit_i_new(rnd_idx)-timestep_size;
% shift other half of 0 durations: delay end time by 5 min
car_trips.anzeit_i_new(car_trips.anzeit_i_new-car_trips.abzeit_i_new==0) = car_trips.anzeit_i_new(car_trips.anzeit_i_new-car_trips.abzeit_i_new==0)+timestep_size;
% create 5 minute indices
car_trips.abzeit_i_ts = round(car_trips.abzeit_i_new/timestep_size);
car_trips.anzeit_i_ts = round(car_trips.anzeit_i_new/timestep_size);
%% Create duration
for i=1:height(car_trips)
    if car_trips.anzeit_i(i) < car_trips.abzeit_i(i)
        car_trips.duration_ts(i) = day_length-car_trips.abzeit_i_ts(i)+car_trips.anzeit_i_ts(i);
    else
        car_trips.duration_ts(i) = car_trips.anzeit_i_ts(i) - car_trips.abzeit_i_ts(i);
    end
end
%% seperate trips by weekday and weekend
trips_wd = car_trips(car_trips.WOTAG <= 5,:);
trips_we = car_trips(car_trips.WOTAG > 5,:);

%create trip duration probability distributions
duration_wd_prob = histcounts(trips_wd{:,"duration_ts"},'BinMethod','integers')/height(trips_wd);
duration_we_prob = histcounts(trips_we{:,"duration_ts"},'BinMethod','integers')/height(trips_we);
%% departure
[N,edges] = histcounts(trips_wd.abzeit_i_ts, 'Normalization','pdf','BinWidth',1);
edges = edges(2:end) - (edges(2)-edges(1))/2;
hold on

plot(edges, movmean(N,30));
[N,edges] = histcounts(trips_we.abzeit_i_ts, 'Normalization','pdf','BinWidth',1);
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, movmean(N,30));
newcolors = [0.75 0 0; 0.266 0.329 0.416; 0.418 0.5117 164/256; 152/256 183/256 225/256];
legend('Weekday', 'Weekend','FontSize',12)
legend boxoff     
ylabel('Probablility')
xlabel('Departure time')
xlim([0,288])
set(gca,'xtick',[1 36 72 108 144 180 216 252 288])
set(gca,'xticklabel',{"00:00" "03:00" "06:00" "09:00" "12:00" "15:00" "18:00" "21:00" "00:00"})
set(gca,'fontname','times', 'FontSize', 12)
set(gcf,'position',[0,0,800,200])
colororder(newcolors)
%% arrival
[N,edges] = histcounts(trips_wd.anzeit_i_ts, 'Normalization','pdf','BinWidth',1);
edges = edges(2:end) - (edges(2)-edges(1))/2;
hold on

plot(edges, movmean(N,30));
[N,edges] = histcounts(trips_we.anzeit_i_ts, 'Normalization','pdf','BinWidth',1);
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, movmean(N,30));
newcolors = [0.75 0 0; 0.266 0.329 0.416; 0.418 0.5117 164/256; 152/256 183/256 225/256];
legend('Weekday', 'Weekend','FontSize',12)
legend boxoff     
ylabel('Probablility')
xlabel('Arrival time')
xlim([0,288])
ylim([0,0.01])
set(gca,'xtick',[1 36 72 108 144 180 216 252 288])
set(gca,'xticklabel',{"00:00" "03:00" "06:00" "09:00" "12:00" "15:00" "18:00" "21:00" "00:00"})
set(gca,'fontname','times', 'FontSize', 12)
set(gcf,'position',[0,0,800,200])
colororder(newcolors)

%Save results

%writematrix([trips_wd.anzeit_i_ts,trips_we.anzeit_i_ts],sprintf('Results/Arrivals_poor.csv'));
writematrix([trips_wd.anzeit_i_ts;trips_we.anzeit_i_ts],sprintf('Results/Arrivals_rich.csv'));
%% Create speed distribution depending on time of day and weekday/weekend
% Differentiate by time of day morning (6-12),afternoon (12-8),night (8-6)
speed_morning_wd=trips_wd(((trips_wd.anzeit_i_new<=720) ),[21,28,29]);
speed_morning_wd=speed_morning_wd((speed_morning_wd.anzeit_i_new> 360),:);

speed_afternoon_wd=trips_wd(((trips_wd.anzeit_i_new<=1200) ),[21,28,29]);
speed_afternoon_wd=speed_afternoon_wd((speed_afternoon_wd.anzeit_i_new> 720),:);

speed_night_wd=vertcat(trips_wd(((trips_wd.anzeit_i_new<=360) ),[21,28,29]),trips_wd(((trips_wd.anzeit_i_new>1200) ),[21,28,29]));

speed_morning_we=trips_we(((trips_we.anzeit_i_new<=720) ),[21,28,29]);
speed_morning_we=speed_morning_we((speed_morning_we.anzeit_i_new> 360),:);

speed_afternoon_we=trips_we(((trips_we.anzeit_i_new<=1200) ),[21,28,29]);
speed_afternoon_we=speed_afternoon_we((speed_afternoon_we.anzeit_i_new> 720),:);

speed_night_we=vertcat(trips_we(((trips_we.anzeit_i_new<=360) ),[21,28,29]),trips_we(((trips_we.anzeit_i_new>1200) ),[21,28,29]));

speed_collected.morning_wd= speed_morning_wd;
speed_collected.afternoon_wd= speed_afternoon_wd;
speed_collected.night_wd= speed_night_wd;
speed_collected.morning_we= speed_morning_we;
speed_collected.afternoon_we = speed_afternoon_we;
speed_collected.night_we= speed_night_we;

hold on  
histogram(speed_morning_wd.SPEED)
histogram(speed_afternoon_wd.SPEED)
histogram(speed_night_wd.SPEED)
hold off

hold on  
histogram(speed_morning_we.SPEED)
histogram(speed_afternoon_we.SPEED)
histogram(speed_night_we.SPEED)
hold off
%% Develop State Chain Based on Trip Data
dataframe = car_trips(:,:);
UUIDs = unique(dataframe.UUID(:).');
no_persons = length(unique(dataframe.UUID));
%create empty table with rows representing timesteps
%and each column representing an individual person
no_timesteps = (60*24*7)/timestep_size;
state_chain = array2table(NaN(no_timesteps,no_persons));
state_chain.Properties.VariableNames = string(UUIDs);
% group by persons 
for i = 1:no_persons
    %sort values according to reported days and reported trip number
    trips = dataframe(dataframe.UUID == UUIDs(i),:);
    trips = sortrows(trips,["BERTAG","WNRT"]);
    %init memory for the prev_status, prev_anzeit, and first reported trips info
    prev_status = -1;
    prev_anzeit_i = -1;
    first_trip = -1;
    %iterate through trips of a person
    for j = 1:height(trips)
        row = trips(j,:);
        weekday_i = weekdayIndex(row.WOTAG);
        abzeit_i = row.abzeit_i_ts + weekday_i;
        anzeit_i = row.anzeit_i_ts + weekday_i;
        %Init status before the first trip on the first reported day
        if prev_status == -1
            first_trip = row;
            % if frist trip on first reported day towards home: Status before this trip = Other
            if row.ZWECK == 7
                state_chain = assignStatus(weekday_i, abzeit_i, weekday_i, string(row.UUID), 3, state_chain);
                first_trip.status = 3;
            %otherwise Status before this trip = Home
            else
                state_chain = assignStatus(weekday_i, abzeit_i, weekday_i, string(row.UUID), 1, state_chain);
                first_trip.status = 1;
            end
        %If no init necessary, assign the status derived from the prev trip until start of this trip
        else  
            state_chain = assignStatus(prev_anzeit_i, abzeit_i, weekday_i, string(row.UUID), prev_status, state_chain);
        end
        % Assign status driving between departure and arrival time depending on means of transport   
        state_chain = assignStatus(abzeit_i, anzeit_i, weekday_i, string(row.UUID), 0, state_chain); 
        %assign memory of the status, anzeit      
        prev_status = returnStatus(row.ZWECK);
        prev_anzeit_i = anzeit_i;
    end
    %Close the loop between last reported trip and first reported trip
    state_chain = assignStatus(prev_anzeit_i, weekdayIndex(first_trip.WOTAG), weekdayIndex(first_trip.WOTAG),string(row.UUID), first_trip.status, state_chain);
end
%% Check if any NaN 
sum(any(isnan(state_chain{:,:})))
%% Extract parking and trip information
% save parking information in parking_matrix
% col 1: location of parking = status
% col 2: parking duration relative for a day (i.e. start counting at 0:00
% and end counting at 23:55
% col 3: day of the week
% col 4: number of trips conducted on this day from this person
parking_idx = 1;
state_chain_array = table2array(state_chain);
% k = persons (each column is person)
for k = 1:width(state_chain_array)
    car_user = false;
    for day = 1:7
        trip_freq = 0;
        prev_stat = -1;
        parking_stat = state_chain_array(1,k);
        parking_duration_idx = 1;
        prev_parking_idx = parking_idx;
        for i = weekdayIndex(day):weekdayIndex(day+1)-1
            %if person starts driving
            if prev_stat ~= state_chain_array(i,k) && state_chain_array(i,k) == 0
                car_user = true;
                prev_stat = state_chain_array(i,k);
                trip_freq = trip_freq + 1;
                parking_matrix(parking_idx,1) = parking_stat;
                parking_matrix(parking_idx,2) = i - parking_duration_idx;
                parking_idx = parking_idx+1;
            %if person changes to a stationary state (i.e. no mean of transport)
            elseif prev_stat ~= state_chain_array(i,k) && state_chain_array(i,k) ~= 0 
                prev_stat = state_chain_array(i,k);
                parking_stat = state_chain_array(i,k);
                parking_duration_idx = i;
            end
        end
        parking_matrix(prev_parking_idx:parking_idx-1,3) = day;
        parking_matrix(prev_parking_idx:parking_idx-1,4) = trip_freq;
    end
end

%% Extract trip frequency for each day of the week and save in frequency_matrix
% col 1: day of the week
% col 2: number of trips on this day of the week
% col 3: flag if the person uses a car
iterator = 1;
frequency_matrix = nan(no_persons*7,1);
for i = 1:no_persons
    trips = W(W.UUID == UUIDs(i),:);
    %iterate through trips of a person
    for day = 0:6
        frequency_matrix(iterator+day,1) = day+1;
        frequency_matrix(iterator+day,2) = height(trips(trips.WOTAG==day+1 & trips.VMDIW==41,1));
        frequency_matrix(iterator+day,3) = any(trips.VMDIW==41);
    end
    iterator = iterator +7;
end
%% Track arrival times of the last trip of the day and fit distributions
% col 1 = anzeit
% col 2 = wotag
car_trips_array = table2array(car_trips(:,["UUID","WOTAG","WNRT","anzeit_i_ts"]));
idx = 1;
for i = 1:no_persons
    %select trips from that person
    trips = car_trips_array(car_trips_array(:,1)==UUIDs(i),:);
    for day = 1:7
        day_trips = trips(trips(:,2)==day,:);
        if ~isempty(day_trips)
            [last_WNRT, last_idx]= max(day_trips(:,3));
            last_trip_anzeit(idx,1) = day_trips(last_idx,4);
            last_trip_anzeit(idx,2) = day_trips(last_idx,2);
            idx = idx+1;
        end
    end
end
wd_arrivals = last_trip_anzeit(last_trip_anzeit(:,2)<=5,:);
we_arrivals = last_trip_anzeit(last_trip_anzeit(:,2)>5,:);
% shift before 7 a.m. last arrivals on tuesday till saturday such that they
% account for weekdays monday to friday, vice versa for sunday and monday
%trips occuring between Tuesday and friday
temp = last_trip_anzeit(last_trip_anzeit(:,2)>=2 | last_trip_anzeit(:,2)<=6,:);
%move early trips to the next day so they are not counted as frist trips of
%the day
temp = temp(temp(:,1)<7*(60/timestep_size),:) + 24*(60/timestep_size);
%start last trip arrivals at 4 pm to fit distribution
wd_arrivals= wd_arrivals(wd_arrivals(:,1)>=16*(60/timestep_size),:);
%then add on late trips (past 12 am) to data set
wd_arrivals(end+1:end+length(temp),:) = temp;
%trips occuring Monday and Sunday
temp = last_trip_anzeit(last_trip_anzeit(:,2)==1 | last_trip_anzeit(:,2)==7,:);
%move early trips to the next day so they are not counted as frist trips of
%the day
temp = temp(temp(:,1)<7*(60/timestep_size),:) + 24*(60/timestep_size);
%start last trip arrivals at 4 pm to fit distribution
we_arrivals= we_arrivals(we_arrivals(:,1)>=16*(60/timestep_size),:);
%then add on late trips (past 12 am) to data set
we_arrivals(end+1:end+length(temp),:) = temp;

%fitdist(wd_arrivals(:,1),"Gamma")
%fitdist(wd_arrivals(:,1),"Lognormal")
%fitdist(wd_arrivals(:,1),"Weibull")
burr_wd = fitdist(wd_arrivals(:,1),"Burr");
burr_we = fitdist(we_arrivals(:,1),"Burr");
%plot(gampdf(1:2000,81.02,13.97))
%plot(lognpdf(1:2000,7.025,0.1096))
%plot(wblpdf(1:2000,1191.67,7.5382))
newcolors = [0.75 0 0; 0.266 0.329 0.416; 0.418 0.5117 164/256; 152/256 183/256 225/256];
hold on
plot(pdf("burr",144:384,burr_we.alpha,burr_we.c,burr_we.k))
plot(pdf("burr",144:384,burr_wd.alpha,burr_wd.c,burr_wd.k))
%plot(pdf("burr",720:1920,1032.87,28.6427,0.32142))
legend('Weekend', 'Weekday','FontSize',12)
legend boxoff   
xlim([1 240])  
ylabel('Probablility')
xlabel('Time')
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'xtick',[1 48 96 144 192 240])
set(gca,'xticklabel',{'12:00' '16:00' '20:00' '24:00' '4:00' '8:00'})
set(gca,'fontname','times', 'FontSize', 12)
set(gcf,'position',[0,0,800,200])
colororder(newcolors)
hold off
%% Seperate data for weekday and weekends
parking_matrix_wd = parking_matrix(parking_matrix(:,3) <= 5,:);
parking_matrix_we = parking_matrix(parking_matrix(:,3) > 5,:);
% Create trip frequency data for cars
car_frequency = frequency_matrix(frequency_matrix(:,3)==1,:);
%% Probability distribution of the states
%calculate the probability of each state at each timestep
probArray = nan(no_timesteps,4);
parfor i = 1:no_timesteps
    probArray(i,:) = histcounts(state_chain_array(i,:),'BinMethod','integers')/no_persons;
end
%% plot state distribution of a weekday
hold off
area(movmean(probArray(864:1152,:),10))
xlim([1 day_length])
ylim([0 1])
set(gca,'fontname','times', 'FontSize', 10)
set(gca,'ytick',[0 1])
set(gca,'yticklabel',{'0' '100'})
%set(gca,'xtick',[1 36 72 108 144 180 216 252 288])
set(gca,'xtick',[1 36:36:day_length])
set(gca,'xticklabel',{'00:00' '03:00' '06:00' '09:00' '12:00' '15:00' '18:00' '21:00' '24:00'})
ylabel('State distribution [%]')
xlabel('Time')
set(gcf,'position',[0,0,800,300])
legend('Driving', 'Home', 'Work', 'Other')
colororder(newcolors)
%% create transition_matrix of shape CxSxS
% Only consider location states: home, work, other
transitionMatrix = nan(no_timesteps, 3, 3);
%determine transition probability at t based on probability to end up at t+1
for i = 2:no_timesteps+1
    %close the loop for the last element
    if i == no_timesteps+1
        prob = probArray(1,:);
    else
        prob = probArray(i,:);
    end
    % iterate through states and normalize the transition probabilities
    % j equals columns therefore the corresponding state coding is j-1
    for j = 1:5
        % if current state home transition towards work/other
        if j == 2
            rel_probs = prob(3)+prob(4);
            transitionMatrix(i-1,1,1) = 0;
            transitionMatrix(i-1,1,2) = prob(3)/rel_probs;
            transitionMatrix(i-1,1,3) = prob(4)/rel_probs;
        %current state work transition towards home/other
        elseif j == 3
            rel_probs = prob(2)+prob(4);
            transitionMatrix(i-1,2,1) = prob(2)/rel_probs;
            transitionMatrix(i-1,2,2) = 0;
            transitionMatrix(i-1,2,3) = prob(4)/rel_probs;
        %current state other transition towards home/work/other
        else
            rel_probs = prob(2)+prob(3)+prob(4);
            transitionMatrix(i-1,3,1) = prob(2)/rel_probs;
            transitionMatrix(i-1,3,2) = prob(3)/rel_probs;
            transitionMatrix(i-1,3,3) = prob(4)/rel_probs;
        end
    end
end
%% Generate trip chain based on markov state transitions and random usage data
% this variable sets the amount of generated chains
iterations = 5000;
rng(6)

samples = nan(no_timesteps,iterations);
test_trips = nan(7,iterations);
tic
for k = 1:iterations
    %simulate for 7 1/2 days to avoid cut off on last day
    chain = nan(no_timesteps+(12*60)/timestep_size,1);
    %init with home
    chain(1) = 1;
    %for each day of the week generate the daily trips
    prev_end_idx = 1;
    for day = 1:7
        % draw the number of trips for that day
        no_trips = randsample(car_frequency(car_frequency(:,1)==day,2),1);
        test_trips(day,k)=no_trips;
        idx = 999999;
        % differ between weekday and weekend
        if day <=5
            param_first_trip = trips_wd(trips_wd.WNRT == 1,:);
            param_parking_matrix = parking_matrix_wd;
            param_duration = duration_wd_prob;
        else
            param_first_trip = trips_we(trips_we.WNRT == 1,:);
            param_parking_matrix = parking_matrix_we;
            param_duration = duration_we_prob;
        end

        % allow next day to be weekend and give some additional time to
        % complete last trip of the day on Friday and Saturday, treat
        % Sunday as weekday
        if day == 5 || day == 6
            endidx = weekdayIndex(day) + random('burr',burr_we.alpha,burr_we.c,burr_we.k);
        else
            endidx = weekdayIndex(day) + random('burr',burr_wd.alpha,burr_wd.c,burr_wd.k);
        end
        iteration = 1;
        %make sure that the trips are completed within a predefined end time

        while 1            
            [chain,idx,tripdurr] = generateDailyTrips(chain,day,param_first_trip,transitionMatrix,param_parking_matrix,no_trips,param_duration,prev_end_idx);
            if (idx <= endidx)                
                prev_end_idx = idx;
                break
            end
            % if endIdx and no_trips are an impossible combination --> redraw no_trips
            iteration = iteration + 1;
            if iteration > 1000
                no_trips = randsample(car_frequency(car_frequency(:,1)==day,2),1);
                test_trips(day,k)=no_trips;
            end
        end
    end
    % fill the gap between last conducted trip and end of week
    chain(idx+1:no_timesteps) = 1;    
    samples(:,k) = chain(1:no_timesteps);
end
toc
sum(any(isnan(samples)))
%% Benchmark average daily driving times of model against data
sumWe = 0;
for k=1:iterations
    for i=(5*24*60)/timestep_size:no_timesteps
        if samples(i,k) == 0
            sumWe = sumWe +timestep_size;
        end
    end
end
sumWd = 0;
for k=1:iterations
    for i=1:-1+(5*24*60)/timestep_size
        if samples(i,k) == 0
            sumWd = sumWd +timestep_size;
        end
    end
end
summ = 0;
for k=1:iterations
    for i=1:no_timesteps
        if samples(i,k) == 0
            summ = summ +timestep_size;
        end
    end
end
disp("MODEL - Average daily driving time total: " + summ/(7*iterations) + ...
    " | weekday: " + sumWd/(5*iterations) + " | weekend: " + sumWe/(2*iterations));
% calculate average trip durations from Dataset
disp("DATA -  Average daily driving time total: " + ...
    (sum(car_trips.duration_ts*timestep_size)/(7*length(unique(car_trips.UUID)))) + ...
    " | weekday: " + sum(trips_wd.duration_ts*timestep_size)/(5*length(unique(car_trips.UUID))) ...
    + " | weekend: " + sum(trips_we.duration_ts*timestep_size)/(2*length(unique(car_trips.UUID))));
disp("MODEL - Average trip frequency: " + sum(sum(test_trips))/(7*iterations)+" | weekday: " + sum(sum(test_trips(1:5,:)))/(5*iterations) + ...
    " | weekend: " + sum(sum(test_trips(6:7,:)))/(2*iterations));
disp("DATA -  Average trip frequency: " + mean(car_frequency(:,2))+" | weekday: " + mean(car_frequency(car_frequency(:,1)<=5,2)) + " | weekend: " + ...
    mean(car_frequency(car_frequency(:,1)>5,2)));
%% Plot the state distributions
% required for visualization (makes sure that each for timestep each state exists)
samples(:,iterations+1) = 0;
samples(:,iterations+2) = 2;
samples(:,iterations+3) = 3;

probSamples = nan(no_timesteps,4);
for i = 1:no_timesteps
    probSamples(i,:) = histcounts(samples(i,:),'BinMethod','integers')/iterations;
end
area(probSamples(864:1152,:))
xlim([1 day_length])
ylim([0 1])
set(gca,'fontname','times', 'FontSize', 12)
set(gca,'ytick',[0 1])
set(gca,'yticklabel',{'0' '100'})
%set(gca,'xtick',[1 36 72 108 144 180 216 252 288])
set(gca,'xtick',[1 36:36:day_length])
set(gca,'xticklabel',{'00:00' '03:00' '06:00' '09:00' '12:00' '15:00' '18:00' '21:00' '24:00'})
ylabel('State distribution [%]')
xlabel('Time')
set(gcf,'position',[0,0,800,300])
%legend('Driving', 'Home', 'Work', 'Other')
colororder(newcolors)
samples = samples(:,1:iterations);
%% Initialize EV DATA


Variable = ["Batterycap"; "Energy Consumption"];
Mini = [32.3;0.177];
Small = [41; 0.203];
Compact = [58;0.193];
Medium = [75;0.209];
SUV = [86.5; 0.258];
Luxury = [83.7; 0.33];
car_data = table(Variable,Mini, Small, Compact, Medium, SUV, Luxury);
speed = speed_collected;
consumption_factor = 1;
charging_power = 11;
charging_eff = 0.9;

%% Preprocess household data
%Share of households with car
%Account for poor/rich HH car proportion

%Test usage
%H=HH_rich;
H16_car = H(H.PKWHH>0,["HHGRO","PKWHH"]);
%% Initialize BESS DATA
BESS_data.capacity = 13.5;
BESS_data.SOC = 50;
BESS_data.charging_efficiency = 0.9;
%% Agglomeration scenarios
% pers per hh: [1; 2; 3; 4 5 or more]
rural.pers_per_hh = [0.35 0.35 0.14 0.12 0.04];
%households per house hh_per_house
rural.hh_per_house = [0.69967127 0.171835405 0.0911751 0.031827008 0.005491218];
rural.hh_per_house_amount = [1 2 5 10 15];
rural.no_nodes = 96;
rural.private_parking_prob = 0.827;
rural.hh_with_car = 0.89;
semi_urb.pers_per_hh = [0.40 0.33 0.12 0.11 0.04];
semi_urb.hh_per_house = [0.552366867 0.131535785 0.160821791 0.118708344 0.036567214];
semi_urb.hh_per_house_amount = [1 2 5 10 15];
semi_urb.no_nodes = 110;
semi_urb.private_parking_prob = 0.825;
semi_urb.hh_with_car = 0.86;
urban.pers_per_hh = [0.54 0.27 0.10 0.07 0.02];
urban.hh_per_house = [0.528080333 0.099991497 0.143230958 0.152232899 0.076464313];
urban.hh_per_house_amount = [1 2 5 10 15];
urban.no_nodes = 112;
urban.private_parking_prob = 0.515;
urban.hh_with_car = 0.77;
%% Set time of year scenarios and charging data

december.consumption_factor = 1.25;
december.household_profiles = householdprofiles_december;
december.name = "December";

march.consumption_factor = 1.15;
march.household_profiles = householdprofiles_march;
march.name = "March";

june.consumption_factor = 1;
june.household_profiles = householdprofiles_june;
june.name = "June";

september.consumption_factor = 1.15;
september.household_profiles = householdprofiles_september;
september.name = "September";


housetype_charger_available.single_house = 0.85;
housetype_charger_available.double_house = 0.72;
housetype_charger_available.multi_house = 0.5;


%% Visualize connection probability
str = '#44546a';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
for current_soc = 0:100
decision(current_soc+1) = min(1,1-1.2/(1+exp(-0.15*(current_soc-60))));
end
plot(decision, 'Color', color)
legend('Connection probability','FontSize',14)
legend boxoff     
xlim([0 100])
ylim([0 1])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'fontname','times', 'FontSize', 12)
xlabel('SOC [%]')

