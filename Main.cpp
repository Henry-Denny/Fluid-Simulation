#include <iostream>
#include <cmath>

#include <SFML/Graphics.hpp>

extern "C"
{
    #include "Fluid.h"
}

#define PI 3.1415296
#define CELL_SIZE 7
#define DIM 256

sf::Color hsv(int hue, float sat, float val);
void Render(sf::RenderWindow *window, FluidSquare *fluid);
void DrawDensity(sf::RenderWindow *window, FluidSquare *fluid);

int main()
{
    sf::RenderWindow window;
    window.create(sf::VideoMode(CELL_SIZE * DIM, CELL_SIZE * DIM), "Fluid Sim");
    
    FluidSquare *fluid = FluidSquareCreate(DIM, 0, 0, 0.2);

    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);
    bool simulating = true;
    bool paused = false;
    while (simulating)
    {
        sf::Event ev;
        while(window.pollEvent(ev))
        {
            if (ev.type == sf::Event::Closed)
            {
                simulating = false;
            }
            if (ev.type == sf::Event::KeyPressed)
            {
                switch (ev.key.code)
                {
                    case (sf::Keyboard::C):
                        // ResetFields(&fluid);
                        break;
                    case (sf::Keyboard::Space):
                        paused = !paused;
                        break;
                    default:
                        break;
                }
            }
        }

        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        sf::Vector2i gridPos {mousePos.x / CELL_SIZE, mousePos.y / CELL_SIZE};

        
        if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
        {
            
            for (int i = -3; i < 3; ++i)
            {
                for (int j = -3; j < 3; ++j)
                {
                    if (gridPos.x + i < 0 || gridPos.y + j < 0 || gridPos.x + i >= fluid->size || gridPos.y + j >= fluid->size) continue;
                    FluidSquareAddDensity(fluid, gridPos.x + i, gridPos.y + j, 50.0f);
                    FluidSquareAddVelocity(fluid, gridPos.x + i, gridPos.y + j, ((float)rand() / RAND_MAX - 0.5f), ((float)rand() / RAND_MAX - 0.5f));
                }
            }
        }
        if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Right))
        {
            for (int i = -3; i < 3; ++i)
            {
                for (int j = -3; j < 3; ++j)
                {
                    sf::Vector2f force = sf::Vector2f(mousePos - prevMousePos);
                    if (force.x != 0.0f || force.y != 0.0f) force = force / std::sqrt(force.x * force.x + force.y * force.y);
                    if (gridPos.x + i < 0 || gridPos.y + j < 0 || gridPos.x + i >= fluid->size || gridPos.y + j >= fluid->size) continue;
                    FluidSquareAddVelocity(fluid, gridPos.x + i, gridPos.y + j, force.x, force.y);
                }
            }
        }
        prevMousePos = mousePos;

        if (!paused)
        {
            FluidSquareAddDensity(fluid, 80, 128, 500.0f);
            FluidSquareAddVelocity(fluid, 80, 128, 0.1f, 0.0f);

            FluidSquareAddDensity(fluid, 176, 128, 500.0f);
            FluidSquareAddVelocity(fluid, 176, 128, -0.1f, 0.0f);

            FluidSquareStep(fluid);
        }

        Render(&window, fluid);
    }

    FluidSquareFree(fluid);
    return 0;
}

void Render(sf::RenderWindow *window, FluidSquare *fluid)
{
    window->clear(sf::Color::Black);

    DrawDensity(window, fluid);

    window->display();
}

void DrawDensity(sf::RenderWindow *window, FluidSquare *fluid)
{
    sf::RectangleShape cell { {CELL_SIZE, CELL_SIZE} };
    int N = fluid->size;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            int d = fluid->density[IX(i,j)];
            cell.setFillColor(hsv(std::min(360, d % 360), 0.8f, std::min(1.0f, d / 255.0f)));
            cell.setPosition(i * CELL_SIZE, j * CELL_SIZE);
            window->draw(cell);
        }
    }
}

sf::Color hsv(int hue, float sat, float val)
{
    hue %= 360;
    while(hue<0) hue += 360;
    
    if(sat<0.f) sat = 0.f;
    if(sat>1.f) sat = 1.f;
    
    if(val<0.f) val = 0.f;
    if(val>1.f) val = 1.f;
    
    int h = hue/60;
    float f = float(hue)/60-h;
    float p = val*(1.f-sat);
    float q = val*(1.f-sat*f);
    float t = val*(1.f-sat*(1-f));
    
    switch(h)
    {
        default:
        case 0:
        case 6: return sf::Color(val*255, t*255, p*255);
        case 1: return sf::Color(q*255, val*255, p*255);
        case 2: return sf::Color(p*255, val*255, t*255);
        case 3: return sf::Color(p*255, q*255, val*255);
        case 4: return sf::Color(t*255, p*255, val*255);
        case 5: return sf::Color(val*255, p*255, q*255);
    }
}