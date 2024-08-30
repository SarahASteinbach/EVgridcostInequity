function [no_persons] = identifyHousehold(str)
%IDENTIFYHOUSEHOLD This function identifiey the number of persons in the
%household according to the households name, Information from
%LoadProfileGenerator
    if contains(str,"CHR09")
        no_persons = 1;
    elseif contains(str,"CHR10")
        no_persons = 1;
    elseif contains(str,"CHR13")
        no_persons = 1;
    elseif contains(str,"CHR14")
        no_persons = 3;
    elseif contains(str,"CHR15")
        no_persons = 6;
    elseif contains(str,"CHR18")
        no_persons = 4;
    elseif contains(str,"CHR19")
        no_persons = 2;
    elseif contains(str,"CHR23")
        no_persons = 1;
    elseif contains(str,"CHR25")
        no_persons = 1;
    elseif contains(str,"CHR26")
        no_persons = 1;
    elseif contains(str,"CHR27")
        no_persons = 4;
    elseif contains(str,"CHR28")
        no_persons = 1;
    elseif contains(str,"CHR29")
        no_persons = 1;
    elseif contains(str,"CHR33")
        no_persons = 2;
    elseif contains(str,"CHR35")
        no_persons = 1;
    elseif contains(str,"CHR37")
        no_persons = 1;
    elseif contains(str,"CHR38")
        no_persons = 1;
    elseif contains(str,"CHR39")
        no_persons = 2;
    elseif contains(str,"CHR40")
        no_persons = 2;
    elseif contains(str,"CHR41")
        no_persons = 5;
    elseif contains(str,"CHR47")
        no_persons = 3;
    elseif contains(str,"CHR48")
        no_persons = 4;
    elseif contains(str,"CHR49")
        no_persons = 3;
    elseif contains(str,"CHR50")
        no_persons = 4;
    elseif contains(str,"CHR51")
        no_persons = 2;
    elseif contains(str,"CHR52")
        no_persons = 1;
    elseif contains(str,"CHR54")
        no_persons = 2;
    elseif contains(str,"CHR58")
        no_persons = 2;
    elseif contains(str,"CHR59")
        no_persons = 5;
    elseif contains(str,"OR01")
        no_persons = 1;
    elseif contains(str,"CHR01")
        no_persons = 2;
    elseif contains(str,"CHR02")
        no_persons = 2;
    elseif contains(str,"CHR03")
        no_persons = 3;
    elseif contains(str,"CHR06")
        no_persons = 1;
    elseif contains(str,"CHR07")
        no_persons = 1;
    elseif contains(str,"CHR08")
        no_persons = 3;
    elseif contains(str,"CHR11")
        no_persons = 1;
    elseif contains(str,"CHR12")
        no_persons = 1;
    elseif contains(str,"CHR16")
        no_persons = 2;
    elseif contains(str,"CHR17")
        no_persons = 2;
    elseif contains(str,"CHR20")
        no_persons = 5;
    elseif contains(str,"CHR21")
        no_persons = 2;
    elseif contains(str,"CHR24")
        no_persons = 1;
    elseif contains(str,"CHR30")
        no_persons = 1;
    elseif contains(str,"CHR31")
        no_persons = 1;
    elseif contains(str,"CHR32")
        no_persons = 2;
    elseif contains(str,"CHR36")
        no_persons = 1;
    elseif contains(str,"CHR42")
        no_persons = 3;
    elseif contains(str,"CHR44")
        no_persons = 4;
    elseif contains(str,"CHR45")
        no_persons = 3;
    elseif contains(str,"CHR53")
        no_persons = 4;
    elseif contains(str,"CHR55")
        no_persons = 2;
    elseif contains(str,"CHR56")
        no_persons = 4;
    elseif contains(str,"CHR57")
        no_persons = 4;
    elseif contains(str,"CHR60")
        no_persons = 3;
    elseif contains(str,"CHR61")
        no_persons = 3;
    elseif contains(str,"CHS01")
        no_persons = 4;
    elseif contains(str,"CHS04")
        no_persons = 2;
    elseif contains(str,"CHS12")
        no_persons = 2;
    elseif contains(str,"CHR05")
        no_persons = 5;
    else
        warning("Number of people not identifyable for household " + str);
        no_persons = str;    
    end
end

